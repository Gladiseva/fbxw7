import glob
import os
import re
import traceback

from valis import registration


CROPS_DIR = "/Users/lollija/phd/fbxw7/crops"
BASE_OUTPUT_DIR = "/Users/lollija/phd/fbxw7/valis_registered_crops"

MARKERS = ["FBXW7", "MYC", "NICD"]
REFERENCE_MARKER = "FBXW7"
CROP_MODE = "overlap"


def numeric_sort_key(value):
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def discover_triplets():
    all_files = glob.glob(os.path.join(CROPS_DIR, "*.ome.tif"))
    triplets = {}

    for path in all_files:
        filename = os.path.basename(path)
        match = re.match(r"^(.+?)_([A-Za-z0-9]+)\.ome\.tif$", filename)
        if not match:
            continue

        sample_id, marker = match.groups()
        marker = marker.upper()
        if marker not in MARKERS:
            continue

        triplets.setdefault(sample_id, {})[marker] = path

    return {
        sample_id: marker_paths
        for sample_id, marker_paths in triplets.items()
        if all(marker in marker_paths for marker in MARKERS)
    }


def register_sample(sample_id, marker_paths):
    reference_path = marker_paths[REFERENCE_MARKER]
    reference_filename = os.path.basename(reference_path)

    slide_paths = [marker_paths[REFERENCE_MARKER]]
    slide_paths.extend(
        marker_paths[marker]
        for marker in MARKERS
        if marker != REFERENCE_MARKER
    )

    results_dst_dir = os.path.join(BASE_OUTPUT_DIR, f"valis_results_sample_{sample_id}")
    registered_slide_dst_dir = os.path.join(results_dst_dir, "registered_slides")
    os.makedirs(registered_slide_dst_dir, exist_ok=True)

    print("")
    print("=" * 70)
    print(f"PROCESSING SAMPLE: {sample_id}")
    print("=" * 70)
    print(f"Reference marker: {REFERENCE_MARKER}")
    print(f"Reference slide: {reference_filename}")
    print("Slides:")
    for path in slide_paths:
        print(f"  - {os.path.basename(path)}")
    print(f"Output folder: {results_dst_dir}")

    registrar = registration.Valis(
        src_dir=CROPS_DIR,
        dst_dir=results_dst_dir,
        img_list=slide_paths,
        reference_img_f=reference_filename,
        align_to_reference=True,
    )

    print(f"[{sample_id}] Running registration...")
    registrar.register()

    print(f"[{sample_id}] Warping and saving registered OME-TIFFs...")
    registrar.warp_and_save_slides(registered_slide_dst_dir, crop=CROP_MODE)
    print(f"[{sample_id}] Done.")


def main():
    print(f"Scanning for OME-TIFF triplets in: {CROPS_DIR}")
    print(f"Using {REFERENCE_MARKER} as the registration reference.")
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    valid_triplets = discover_triplets()
    sample_ids = sorted(valid_triplets, key=numeric_sort_key)

    print(f"Found {len(sample_ids)} complete triplets:")
    print(sample_ids)

    failed_samples = []

    try:
        for sample_id in sample_ids:
            try:
                register_sample(sample_id, valid_triplets[sample_id])
            except Exception as exc:
                print(f"Error processing sample {sample_id}: {exc}")
                traceback.print_exc()
                failed_samples.append(sample_id)
                continue
    finally:
        print("")
        print("=" * 70)
        print("Terminating Java Virtual Machine...")
        registration.kill_jvm()
        print("JVM closed.")

    print("")
    print("=" * 70)
    if failed_samples:
        print(f"Finished with failures: {failed_samples}")
    else:
        print("All crop triplets registered successfully.")
    print("=" * 70)


if __name__ == "__main__":
    main()
