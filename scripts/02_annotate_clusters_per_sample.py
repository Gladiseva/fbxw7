import gc
import glob
import os

import anndata as ad
import napari
import numpy as np
import pandas as pd
import tifffile
from magicgui.widgets import ComboBox, Container, Label, PushButton

from config_utils import get_markers, load_config, parse_config_arg, resolve_path


args = parse_config_arg("Annotate per-sample Leiden clusters in Napari.")
config = load_config(args.config)

RESULTS_DIR = resolve_path(config, "results_per_sample_dir")
CROPS_DIR = resolve_path(config, "crops_dir")
PATCH_SIZE = config["analysis"]["patch_size"]

MARKERS = get_markers(config)
ANNOTATIONS = list(config["annotation_ui"]["labels"])
CLUSTER_COLORS = list(config["annotation_ui"]["cluster_colors"])


def numeric_sort_key(value):
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def standardize_image(img_array):
    img_array = np.squeeze(img_array)

    if img_array.ndim == 3 and img_array.shape[0] in [1, 3, 4]:
        img_array = np.transpose(img_array, (1, 2, 0))
    if img_array.ndim == 2:
        img_array = np.stack((img_array,) * 3, axis=-1)
    if img_array.ndim == 3 and img_array.shape[-1] == 4:
        img_array = img_array[..., :3]

    return img_array


def patch_rectangles(obs_df):
    rectangles = []
    for _, row in obs_df.iterrows():
        y = int(row["y"])
        x = int(row["x"])
        rectangles.append(
            np.array(
                [
                    [y, x],
                    [y, x + PATCH_SIZE],
                    [y + PATCH_SIZE, x + PATCH_SIZE],
                    [y + PATCH_SIZE, x],
                ]
            )
        )
    return rectangles


class PerSampleClusterAnnotator:
    def __init__(self, viewer):
        self.viewer = viewer

        self.marker = None
        self.samples = []
        self.curr_sample_idx = 0
        self.curr_cluster_idx = 0
        self.clusters = []

        self.adata = None
        self.adata_path = ""
        self.dirty = False

        self.marker_combo = ComboBox(choices=MARKERS, label="1. Marker:")
        self.load_btn = PushButton(text="Load Marker")

        self.status_label = Label(value="Status: Ready.")
        self.sample_label = Label(value="Sample: -- (0/0)")
        self.cluster_label = Label(value="Cluster: -- (0/0)")

        self.prev_sample_btn = PushButton(text="Prev Sample")
        self.next_sample_btn = PushButton(text="Next Sample")
        self.prev_cluster_btn = PushButton(text="Prev Cluster")
        self.next_cluster_btn = PushButton(text="Next Cluster")

        self.annot_combo = ComboBox(choices=ANNOTATIONS, label="2. Annotation:")
        self.save_annot_btn = PushButton(text="Save Label & Next Cluster")
        self.save_sample_btn = PushButton(text="Save Current Sample")

        self.widget = Container(
            widgets=[
                self.marker_combo,
                self.load_btn,
                Label(value="-" * 40),
                self.status_label,
                self.sample_label,
                self.cluster_label,
                Container(
                    widgets=[self.prev_sample_btn, self.next_sample_btn],
                    layout="horizontal",
                ),
                Container(
                    widgets=[self.prev_cluster_btn, self.next_cluster_btn],
                    layout="horizontal",
                ),
                Label(value="-" * 40),
                self.annot_combo,
                self.save_annot_btn,
                self.save_sample_btn,
            ]
        )

        self.load_btn.clicked.connect(self.load_marker)
        self.prev_sample_btn.clicked.connect(self.prev_sample)
        self.next_sample_btn.clicked.connect(self.next_sample)
        self.prev_cluster_btn.clicked.connect(self.prev_cluster)
        self.next_cluster_btn.clicked.connect(self.next_cluster)
        self.save_annot_btn.clicked.connect(self.save_annotation)
        self.save_sample_btn.clicked.connect(self.save_to_disk)

    def load_marker(self):
        self.save_if_dirty()

        self.marker = self.marker_combo.value
        paths = glob.glob(os.path.join(RESULTS_DIR, f"adata_*_{self.marker}.h5ad"))
        self.samples = sorted(
            [
                os.path.basename(path)
                .replace("adata_", "")
                .replace(f"_{self.marker}.h5ad", "")
                for path in paths
            ],
            key=numeric_sort_key,
        )

        if not self.samples:
            self.status_label.value = f"No per-sample AnnData files for {self.marker}."
            self.viewer.layers.clear()
            return

        self.curr_sample_idx = 0
        self.curr_cluster_idx = 0
        self.load_sample()

    def load_sample(self):
        if not self.marker or not self.samples:
            return

        self.save_if_dirty()

        sample = self.samples[self.curr_sample_idx]
        self.adata_path = os.path.join(RESULTS_DIR, f"adata_{sample}_{self.marker}.h5ad")
        self.status_label.value = f"Loading sample {sample} / {self.marker}..."

        if not os.path.exists(self.adata_path):
            self.status_label.value = f"Missing AnnData: {self.adata_path}"
            return

        self.adata = ad.read_h5ad(self.adata_path)
        self.dirty = False

        required_cols = {"leiden", "x", "y"}
        missing = required_cols.difference(self.adata.obs.columns)
        if missing:
            self.status_label.value = f"Missing columns in AnnData: {sorted(missing)}"
            return

        self.ensure_annotation_column()
        self.clusters = sorted(
            self.adata.obs["leiden"].astype(str).unique(),
            key=numeric_sort_key,
        )
        self.curr_cluster_idx = min(self.curr_cluster_idx, max(len(self.clusters) - 1, 0))

        self.render_sample()

    def ensure_annotation_column(self):
        if "annotation" not in self.adata.obs.columns:
            values = ["unannotated"] * self.adata.n_obs
        else:
            values = self.adata.obs["annotation"].astype(str).fillna("unannotated").tolist()
            values = [value if value in ANNOTATIONS else "unannotated" for value in values]

        self.adata.obs["annotation"] = pd.Categorical(values, categories=ANNOTATIONS)

    def render_sample(self):
        if self.adata is None:
            return

        sample = self.samples[self.curr_sample_idx]
        slide_path = os.path.join(CROPS_DIR, f"{sample}_{self.marker}.ome.tif")
        if not os.path.exists(slide_path):
            self.status_label.value = f"Missing crop image: {slide_path}"
            return

        self.status_label.value = "Rendering sample..."
        self.viewer.layers.clear()

        try:
            with tifffile.TiffFile(slide_path) as tif:
                img_array = standardize_image(tif.series[0].asarray())
            self.viewer.add_image(img_array, name=f"{sample}_{self.marker}")
            del img_array
            gc.collect()
        except Exception as exc:
            self.status_label.value = f"Error loading image: {exc}"
            return

        all_rectangles = []
        edge_colors = []
        leiden_as_str = self.adata.obs["leiden"].astype(str)
        for cluster_idx, cluster in enumerate(self.clusters):
            mask = leiden_as_str == cluster
            rectangles = patch_rectangles(self.adata.obs.loc[mask])
            all_rectangles.extend(rectangles)
            edge_colors.extend([CLUSTER_COLORS[cluster_idx % len(CLUSTER_COLORS)]] * len(rectangles))

        if all_rectangles:
            self.viewer.add_shapes(
                all_rectangles,
                shape_type="rectangle",
                edge_color=edge_colors,
                face_color="transparent",
                edge_width=2,
                name="All sample clusters",
            )

        self.update_highlight()

    def update_highlight(self):
        if self.adata is None or not self.clusters:
            return

        for layer in list(self.viewer.layers):
            if layer.name == "Current cluster":
                self.viewer.layers.remove(layer)

        sample = self.samples[self.curr_sample_idx]
        cluster = self.clusters[self.curr_cluster_idx]
        leiden_as_str = self.adata.obs["leiden"].astype(str)
        mask = leiden_as_str == cluster
        subset = self.adata.obs.loc[mask]
        rectangles = patch_rectangles(subset)

        if rectangles:
            self.viewer.add_shapes(
                rectangles,
                shape_type="rectangle",
                edge_color="yellow",
                face_color="transparent",
                edge_width=7,
                name="Current cluster",
            )

        saved_values = subset["annotation"].astype(str).unique().tolist()
        saved_values = [value for value in saved_values if value in ANNOTATIONS]
        if len(saved_values) == 1:
            self.annot_combo.value = saved_values[0]
        elif "unannotated" in ANNOTATIONS:
            self.annot_combo.value = "unannotated"

        self.sample_label.value = (
            f"Sample: {sample} ({self.curr_sample_idx + 1} / {len(self.samples)})"
        )
        self.cluster_label.value = (
            f"Cluster: {cluster} ({self.curr_cluster_idx + 1} / {len(self.clusters)})"
        )
        self.status_label.value = (
            f"Ready. Current cluster has {len(subset)} patches."
        )

    def prev_sample(self):
        if not self.samples or self.curr_sample_idx == 0:
            return
        self.curr_sample_idx -= 1
        self.curr_cluster_idx = 0
        self.load_sample()

    def next_sample(self):
        if not self.samples or self.curr_sample_idx >= len(self.samples) - 1:
            self.status_label.value = "Last sample for this marker."
            return
        self.curr_sample_idx += 1
        self.curr_cluster_idx = 0
        self.load_sample()

    def prev_cluster(self):
        if not self.clusters or self.curr_cluster_idx == 0:
            return
        self.curr_cluster_idx -= 1
        self.update_highlight()

    def next_cluster(self):
        if not self.clusters:
            return
        if self.curr_cluster_idx < len(self.clusters) - 1:
            self.curr_cluster_idx += 1
            self.update_highlight()
        else:
            self.status_label.value = "Last cluster in this sample."

    def save_annotation(self):
        if self.adata is None or not self.clusters:
            return

        cluster = self.clusters[self.curr_cluster_idx]
        annot_val = self.annot_combo.value
        mask = self.adata.obs["leiden"].astype(str) == cluster
        self.adata.obs.loc[mask, "annotation"] = annot_val
        self.dirty = True
        self.save_to_disk()

        if self.curr_cluster_idx < len(self.clusters) - 1:
            self.curr_cluster_idx += 1
            self.update_highlight()
            return

        if self.curr_sample_idx < len(self.samples) - 1:
            self.curr_sample_idx += 1
            self.curr_cluster_idx = 0
            self.load_sample()
            return

        self.status_label.value = "All samples completed for this marker."

    def save_if_dirty(self):
        if self.dirty and self.adata is not None:
            self.save_to_disk()

    def save_to_disk(self):
        if self.adata is None or not self.adata_path:
            return

        self.status_label.value = "Saving current sample..."
        self.adata.write_h5ad(self.adata_path)
        self.dirty = False
        self.status_label.value = f"Saved: {os.path.basename(self.adata_path)}"


if __name__ == "__main__":
    viewer = napari.Viewer(title="Per-Sample Cluster Annotator")
    annotator = PerSampleClusterAnnotator(viewer)
    viewer.window.add_dock_widget(
        annotator.widget,
        name="Per-Sample Annotation Toolkit",
        area="right",
    )
    napari.run()
