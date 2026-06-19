import glob
import os
import re
from collections import defaultdict

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # <-- NEW: Imported pandas to handle CSV saving
import squidpy as sq

IN_DIR = "/Users/lollija/phd/fbxw7/results_per_sample_valis_coords"
OUT_STATS_DIR = "/Users/lollija/phd/fbxw7/spatial_stats_results_marker_status"

# Patch size and spatial graph parameters
PATCH_SIZE = 224
RADIUS = PATCH_SIZE * 2  

def parse_adata_path(path):
    filename = os.path.basename(path)
    match = re.match(r"^adata_(.+?)_([A-Za-z0-9]+)\.h5ad$", filename)
    if not match:
        return None
    return match.groups()

def prepare_spatial_coords(adata):
    if "valis_registered_centroid_x" not in adata.obs:
        raise ValueError("VALIS registered coordinates not found in adata.obs.")
    
    coords = np.column_stack([
        adata.obs["valis_registered_centroid_x"].to_numpy(),
        adata.obs["valis_registered_centroid_y"].to_numpy()
    ])
    adata.obsm["spatial"] = coords
    return adata

def main():
    os.makedirs(OUT_STATS_DIR, exist_ok=True)

    adata_paths = glob.glob(os.path.join(IN_DIR, "adata_*.h5ad"))
    files_by_sample = defaultdict(list)
    
    for path in adata_paths:
        parsed = parse_adata_path(path)
        if parsed:
            sample_id, marker = parsed
            files_by_sample[sample_id].append((marker, path))

    all_sample_adatas = []

    print(f"Found {len(files_by_sample)} samples. Beginning marker-status spatial analysis...")

    for sample_id, marker_files in files_by_sample.items():
        print(f"\n--- Processing Sample: {sample_id} ---")
        
        sample_marker_adatas = []
        for marker, path in marker_files:
            adata = ad.read_h5ad(path)
            
            # FIX WARNING: Make obs_names globally unique by appending both sample_id and marker
            adata.obs_names = adata.obs_names + f"_{sample_id}_{marker}"
            adata.obs['marker'] = marker
            adata.obs['sample_id'] = sample_id
            
            try:
                adata = prepare_spatial_coords(adata)
                sample_marker_adatas.append(adata)
            except ValueError as e:
                print(f"Skipping {marker} in {sample_id}: {e}")
                continue
        
        if not sample_marker_adatas:
            continue
            
        adata_sample = ad.concat(sample_marker_adatas, join="outer")
        
        # --- FILTER OUT 'NOT TUMOR' PATCHES ---
        if "annotation" in adata_sample.obs.columns:
            adata_sample = adata_sample[adata_sample.obs["annotation"] != "not tumor"].copy()
            adata_sample.obs["marker_status"] = adata_sample.obs["marker"].astype(str) + "_" + adata_sample.obs["annotation"].astype(str)
        else:
            print(f"Warning: 'annotation' missing for {sample_id}. Using marker name only.")
            adata_sample.obs["marker_status"] = adata_sample.obs["marker"].astype(str)

        if len(adata_sample) == 0:
            print(f"Sample {sample_id} is empty after filtering out 'not tumor' patches. Skipping.")
            continue

        # Clean categorical data and set types properly
        adata_sample.obs["marker_status"] = adata_sample.obs["marker_status"].astype("category")
        if hasattr(adata_sample.obs["marker_status"].cat, "remove_unused_categories"):
            adata_sample.obs["marker_status"] = adata_sample.obs["marker_status"].cat.remove_unused_categories()
            
        adata_sample.obs["sample_id"] = adata_sample.obs["sample_id"].astype("category")

        # --- PER SAMPLE SQUIDPY ANALYSIS ---
        print(f"Building spatial graph (radius={RADIUS}) for {sample_id}...")
        sq.gr.spatial_neighbors(
            adata_sample, 
            radius=RADIUS, 
            coord_type="generic", 
            spatial_key="spatial"
        )
        
        print("Calculating neighborhood enrichment for marker states...")
        sq.gr.nhood_enrichment(adata_sample, cluster_key="marker_status")
        
        # Replace NaNs/Infs
        zscores = adata_sample.uns["marker_status_nhood_enrichment"]["zscore"]
        cleaned_zscores = np.nan_to_num(zscores, nan=0.0, posinf=0.0, neginf=0.0)
        adata_sample.uns["marker_status_nhood_enrichment"]["zscore"] = cleaned_zscores
        
        # --- NEW: SAVE PER-SAMPLE Z-SCORES TO CSV ---
        categories = adata_sample.obs["marker_status"].cat.categories
        df_zscores = pd.DataFrame(cleaned_zscores, index=categories, columns=categories)
        csv_path = os.path.join(OUT_STATS_DIR, f"marker_status_zscores_sample_{sample_id}.csv")
        df_zscores.to_csv(csv_path)
        print(f"Saved Z-score CSV for {sample_id}.")

        # Plot
        plot_path = os.path.join(OUT_STATS_DIR, f"marker_status_sample_{sample_id}.png")
        sq.pl.nhood_enrichment(
            adata_sample, 
            cluster_key="marker_status", 
            method="average", 
            cmap="inferno",
            title=f"Marker Status Interactions - Sample {sample_id}",
            figsize=(10, 10)
        )
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()
        
        all_sample_adatas.append(adata_sample)

    # --- GLOBAL STATISTICS ---
    if all_sample_adatas:
        print("\n--- Processing GLOBAL Statistics ---")
        adata_global = ad.concat(all_sample_adatas, join="outer")
        
        adata_global.obs["marker_status"] = adata_global.obs["marker_status"].astype("category")
        if hasattr(adata_global.obs["marker_status"].cat, "remove_unused_categories"):
            adata_global.obs["marker_status"] = adata_global.obs["marker_status"].cat.remove_unused_categories()
            
        adata_global.obs["sample_id"] = adata_global.obs["sample_id"].astype("category")

        print("Building global spatial graph across all samples...")
        sq.gr.spatial_neighbors(
            adata_global, 
            radius=RADIUS, 
            coord_type="generic", 
            spatial_key="spatial",
            library_key="sample_id" 
        )
        
        print("Calculating global neighborhood enrichment...")
        sq.gr.nhood_enrichment(adata_global, cluster_key="marker_status")
        
        # Replace NaNs/Infs
        zscores_global = adata_global.uns["marker_status_nhood_enrichment"]["zscore"]
        cleaned_zscores_global = np.nan_to_num(zscores_global, nan=0.0, posinf=0.0, neginf=0.0)
        adata_global.uns["marker_status_nhood_enrichment"]["zscore"] = cleaned_zscores_global
        
        # --- NEW: SAVE GLOBAL Z-SCORES TO CSV ---
        categories_global = adata_global.obs["marker_status"].cat.categories
        df_zscores_global = pd.DataFrame(cleaned_zscores_global, index=categories_global, columns=categories_global)
        global_csv_path = os.path.join(OUT_STATS_DIR, "marker_status_zscores_GLOBAL.csv")
        df_zscores_global.to_csv(global_csv_path)
        print("Saved Global Z-score CSV.")

        # Plot
        plot_path = os.path.join(OUT_STATS_DIR, "marker_status_GLOBAL.png")
        sq.pl.nhood_enrichment(
            adata_global, 
            cluster_key="marker_status", 
            method="average", 
            cmap="inferno",
            title="Global Marker Status Interactions (All Samples)",
            figsize=(12, 12)
        )
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()
        
        global_out_path = os.path.join(OUT_STATS_DIR, "adata_global_marker_status.h5ad")
        adata_global.write_h5ad(global_out_path)
        print(f"Global AnnData saved to: {global_out_path}")

    print("\nDone! Check your plots and CSVs in:", OUT_STATS_DIR)

if __name__ == "__main__":
    main()