import glob
import os
import re
from collections import defaultdict

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
from valis import registration

from config_utils import get_markers, load_config, parse_config_arg, resolve_path


args = parse_config_arg("Warp AnnData patch centroids into VALIS-registered coordinates.")
config = load_config(args.config)

ADATA_DIR = resolve_path(config, "results_per_sample_dir")
CROPS_DIR = resolve_path(config, "crops_dir")
VALIS_DIR = resolve_path(config, "valis_registered_crops_dir")
OUT_DIR = resolve_path(config, "results_per_sample_valis_coords_dir")
PLOTS_DIR = resolve_path(config, "results_per_sample_valis_plots_dir")

MARKERS = get_markers(config)
PATCH_SIZE = config["analysis"]["patch_size"]

# Consistent plotting colors for your specific markers
MARKER_COLORS = dict(config["valis_warp"]["marker_colors"])

# Match the setting used in register_crops_valis.py when saving registered slides.
CROP = config["valis_warp"]["crop"]
USE_NON_RIGID = config["valis_warp"]["use_non_rigid"]


def numeric_sort_key(value):
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def parse_adata_path(path):
    filename = os.path.basename(path)
    match = re.match(r"^adata_(.+?)_([A-Za-z0-9]+)\.h5ad$", filename)
    if not match:
        return None

    sample_id, marker = match.groups()
    marker = marker.upper()
    if marker not in MARKERS:
        return None

    return sample_id, marker


def get_registrar_path(sample_id):
    return os.path.join(
        VALIS_DIR,
        f"valis_results_sample_{sample_id}",
        "crops",
        "data",
        "crops_registrar.pickle",
    )


def get_slide(registrar, sample_id, marker):
    candidates = [
        f"{sample_id}_{marker}.ome.tif",
        f"{sample_id}_{marker}",
        os.path.join(CROPS_DIR, f"{sample_id}_{marker}.ome.tif"),
    ]

    for candidate in candidates:
        try:
            return registrar.get_slide(candidate)
        except Exception:
            continue

    available = sorted(getattr(registrar, "slide_dict", {}).keys())
    raise KeyError(
        f"Could not find VALIS slide for sample={sample_id}, marker={marker}. "
        f"Available slide keys: {available}"
    )


def add_registered_centroids(adata, slide_obj):
    x = adata.obs["x"].astype(float).to_numpy()
    y = adata.obs["y"].astype(float).to_numpy()

    centroid_x = x + PATCH_SIZE / 2
    centroid_y = y + PATCH_SIZE / 2
    xy = np.column_stack([centroid_x, centroid_y])

    warped_xy = slide_obj.warp_xy(
        xy,
        slide_level=0,
        pt_level=0,
        non_rigid=USE_NON_RIGID,
        crop=CROP,
    )

    adata.obs["centroid_x"] = centroid_x
    adata.obs["centroid_y"] = centroid_y
    adata.obs["valis_registered_centroid_x"] = warped_xy[:, 0]
    adata.obs["valis_registered_centroid_y"] = warped_xy[:, 1]
    adata.obs["valis_registration_crop"] = CROP
    adata.obs["valis_registration_non_rigid"] = USE_NON_RIGID

    return adata, xy, warped_xy


def plot_registration_results(sample_id, original_coords, registered_coords, out_dir):
    """Generates and saves a side-by-side plot comparing original vs registered coordinates."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Plot original coordinates
    for marker, coords in original_coords.items():
        color = MARKER_COLORS.get(marker, "tab:gray")
        axes[0].scatter(coords[:, 0], coords[:, 1], s=4, label=marker, alpha=0.5, color=color)
    
    axes[0].set_title(f"Sample {sample_id} - Original Unregistered Centroids")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    axes[0].invert_yaxis()  # Image coordinates typically have (0,0) at the top-left
    axes[0].axis('equal')
    axes[0].legend(loc="upper right")

    # Plot VALIS registered coordinates
    for marker, coords in registered_coords.items():
        color = MARKER_COLORS.get(marker, "tab:gray")
        axes[1].scatter(coords[:, 0], coords[:, 1], s=4, label=marker, alpha=0.5, color=color)
        
    axes[1].set_title(f"Sample {sample_id} - VALIS Registered Centroids")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    axes[1].invert_yaxis()
    axes[1].axis('equal')
    axes[1].legend(loc="upper right")

    plt.tight_layout()
    plot_path = os.path.join(out_dir, f"registration_plot_sample_{sample_id}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)

    adata_paths = sorted(
        glob.glob(os.path.join(ADATA_DIR, "adata_*.h5ad")),
        key=lambda p: parse_adata_path(p) or ("", ""),
    )

    # Group files by sample_id to allow per-sample plotting
    files_by_sample = defaultdict(list)
    for path in adata_paths:
        parsed = parse_adata_path(path)
        if parsed:
            sample_id, marker = parsed
            files_by_sample[sample_id].append((marker, path))

    registrar_cache = {}
    skipped = []
    written = []
    plots_written = []

    try:
        for sample_id, marker_files in files_by_sample.items():
            registrar_path = get_registrar_path(sample_id)
            if not os.path.exists(registrar_path):
                skipped.append((sample_id, "ALL", "missing VALIS registrar"))
                continue

            if sample_id not in registrar_cache:
                print(f"Loading VALIS registrar for sample {sample_id}...")
                registrar_cache[sample_id] = registration.load_registrar(registrar_path)

            registrar = registrar_cache[sample_id]
            
            sample_orig_coords = {}
            sample_reg_coords = {}

            for marker, adata_path in marker_files:
                print(f"  Warping centroids for {sample_id}_{marker}...")
                
                try:
                    slide_obj = get_slide(registrar, sample_id, marker)
                except KeyError as e:
                    skipped.append((sample_id, marker, str(e)))
                    continue

                adata = ad.read_h5ad(adata_path)
                adata, xy, warped_xy = add_registered_centroids(adata, slide_obj)

                # Save coordinates for the per-sample plot
                sample_orig_coords[marker] = xy
                sample_reg_coords[marker] = warped_xy

                # Save the AnnData file
                out_path = os.path.join(OUT_DIR, os.path.basename(adata_path))
                adata.write_h5ad(out_path)
                written.append(out_path)

                del adata
            
            # Generate the plot once per sample if we processed any markers successfully
            if sample_orig_coords and sample_reg_coords:
                print(f"  Generating overlay plot for sample {sample_id}...")
                plot_registration_results(sample_id, sample_orig_coords, sample_reg_coords, PLOTS_DIR)
                plots_written.append(sample_id)

    finally:
        registration.kill_jvm()

    print("\n--- Summary ---")
    print(f"Wrote {len(written)} AnnData files with VALIS-registered centroids:")
    for path in written:
        print(f"  - {path}")

    print(f"\nCreated overlay plots for {len(plots_written)} samples in:")
    print(f"  {PLOTS_DIR}")

    if skipped:
        print(f"\nSkipped {len(skipped)} items:")
        for sample_id, marker, reason in skipped:
            print(f"  - {sample_id}_{marker}: {reason}")


if __name__ == "__main__":
    main()
