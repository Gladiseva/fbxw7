# FBXW7 Manuscript Analysis Scripts

This repository contains the analysis scripts accompanying the manuscript. The scripts are kept as a numbered pipeline.

Large input data and generated results are not tracked in git. Place them in the folders shown below, or edit `config/default.yaml` to point to another location.

## Repository Layout

```text
annotation/                         GeoJSON tissue annotations, one per sample-marker image
crops/                              OME-TIFF crop images
results_per_sample/                 Feature, DAB, PCA, Leiden, and manual annotation AnnData files
valis_registered_crops/             VALIS registration outputs
results_per_sample_valis_coords/    AnnData files with registered centroids
results_per_sample_valis_plots/     Registration QC plots
spatial_stats_results_marker_status/Spatial enrichment plots, matrices, and cohort statistics
scripts/                            Numbered analysis scripts
config/default.yaml                 Default paths and analysis parameters
```

Expected input filenames:

```text
crops/{sample_id}_{marker}.ome.tif
annotation/{sample_id}_{marker}.geojson
```

The default markers are `FBXW7`, `MYC`, and `NICD`.

## Environment

Create the environment with:

```bash
conda env create -f environment.yml
conda activate fbxw7
```

## Configuration

All scripts accept the same optional config argument:

```bash
python scripts/01_extract_fm_features.py --config config/default.yaml
```

Relative paths in `config/default.yaml` are resolved against `project_root`. The default `project_root: ..` means paths are relative to the repository root because the config file lives in `config/`.

## Pipeline

Run the scripts in this order:

```bash
python scripts/01_extract_fm_features.py --config config/default.yaml
python scripts/02_annotate_clusters_per_sample.py --config config/default.yaml
python scripts/04_register_crops_valis.py --config config/default.yaml
python scripts/05_warp_adata_centroids_with_valis.py --config config/default.yaml
python scripts/06_calculate_spatial_stats.py --config config/default.yaml
python scripts/07_stats_significance.py --config config/default.yaml
```

`03_exploration.ipynb` is retained as an exploratory notebook.

## Notes

- `01_extract_fm_features.py` skips existing AnnData outputs, so it can be resumed.
- `02_annotate_clusters_per_sample.py` opens a Napari interface and writes annotations back into the per-sample AnnData files.
- Generated data and figures are ignored by git; see `.gitignore`.
