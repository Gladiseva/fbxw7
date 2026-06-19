import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["PYTORCH_MPS_HIGH_WATERMARK_RATIO"] = "0.0"

import gc
import json
import glob
import cv2
import numpy as np
import torch
import timm
import tifffile
import anndata as ad
import scanpy as sc
from PIL import Image
from torchvision import transforms
from tqdm.auto import tqdm
from skimage.color import rgb2hed

# --- Paths ---
CROPS_DIR = "/Users/lollija/phd/fbxw7/crops"
ANN_DIR   = "/Users/lollija/phd/fbxw7/annotation"
OUTPUT_DIR = "/Users/lollija/phd/fbxw7/results_per_sample"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Hyperparameters ---
MARKERS           = ["FBXW7", "MYC", "NICD"]
PATCH_SIZE        = 224
BATCH_SIZE        = 32
BG_THRESHOLD      = 250
MASK_OVERLAP_FRAC = 0.50  # At least 50% of the patch must be inside the GeoJSON mask

# Get all unique sample IDs by parsing the file names in the crops directory
crop_files = [f for f in os.listdir(CROPS_DIR) if f.endswith(".ome.tif")]
sample_ids = sorted(list(set([f.split('_')[0] for f in crop_files])))
print(f"Found {len(sample_ids)} unique samples: {sample_ids}")

device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print(f"Using device: {device}")

model = timm.create_model('hf_hub:JWonderLand/StainNet-Base', pretrained=True)
model = model.to(device)
model.eval()

preprocess = transforms.Compose([
    transforms.Resize(224, interpolation=transforms.InterpolationMode.BICUBIC),
    transforms.ToTensor(),
    transforms.Normalize(mean=(0.5, 0.5, 0.5), std=(0.5, 0.5, 0.5)),
])

def geojson_to_mask(geojson_path, height, width):
    """Rasterizes a GeoJSON into a binary numpy mask (1=tissue, 0=background)"""
    with open(geojson_path, 'r') as f:
        data = json.load(f)
        
    mask = np.zeros((height, width), dtype=np.uint8)
    
    for feat in data.get('features', []):
        geom = feat.get('geometry', {})
        geom_type = geom.get('type', '')
        coords = geom.get('coordinates', [])
        
        if not coords: continue
            
        if geom_type == 'Polygon':
            pts = np.array(coords[0], np.int32)
            cv2.fillPoly(mask, [pts], 1)
        elif geom_type == 'MultiPolygon':
            for poly in coords:
                pts = np.array(poly[0], np.int32)
                cv2.fillPoly(mask, [pts], 1)
                
    return mask

for marker in MARKERS:
    print(f"\n{'='*50}\n▶ Processing Marker: {marker}\n{'='*50}")
    
    for sample_id in sample_ids:
        slide_name   = f"{sample_id}_{marker}"
        slide_path   = os.path.join(CROPS_DIR, f"{slide_name}.ome.tif")
        geojson_path = os.path.join(ANN_DIR, f"{slide_name}.geojson")
        adata_out    = os.path.join(OUTPUT_DIR, f"adata_{slide_name}.h5ad")
        
        # Check if output already exists to allow resuming
        if os.path.exists(adata_out):
            print(f"⏭  Skipping {slide_name} (AnnData exists)")
            continue
            
        if not os.path.exists(slide_path) or not os.path.exists(geojson_path):
            print(f"⚠️  Missing image or annotation for {slide_name}, skipping.")
            continue
            
        print(f"\n🔬 Sample: {slide_name}")
        
        # 1. Open TIF & generate Mask
        tif = tifffile.TiffFile(slide_path)
        series = tif.series[0]
        axes = series.axes
        img_array = series.asarray()
        
        y_idx = axes.find('Y')
        x_idx = axes.find('X')
        slide_h = img_array.shape[y_idx]
        slide_w = img_array.shape[x_idx]
        
        mask = geojson_to_mask(geojson_path, slide_h, slide_w)
        
        # 2. Extract Patches & Features
        coords_grid = [
            (left, top)
            for top  in range(0, slide_h - PATCH_SIZE + 1, PATCH_SIZE)
            for left in range(0, slide_w  - PATCH_SIZE + 1, PATCH_SIZE)
        ]
        
        features_list, dab_list, valid_coords = [], [], []
        buf = []

        def flush_batch(buf_local):
            tensors = torch.stack([b[0] for b in buf_local])
            with torch.autocast(device_type="mps", dtype=torch.float16):
                with torch.inference_mode():
                    feats = model(tensors.to(device)).cpu()
            return feats, [b[1] for b in buf_local], [b[2] for b in buf_local]

        for (left, top) in tqdm(coords_grid, desc="Extracting", leave=False):
            # Mask check: Ensure patch is sufficiently inside the annotation
            patch_mask = mask[top:top+PATCH_SIZE, left:left+PATCH_SIZE]
            if patch_mask.mean() < MASK_OVERLAP_FRAC:
                continue
                
            # Grab patch
            fetch_slices = [slice(None)] * img_array.ndim
            fetch_slices[y_idx] = slice(top, top + PATCH_SIZE)
            fetch_slices[x_idx] = slice(left, left + PATCH_SIZE)
            
            patch_np = img_array[tuple(fetch_slices)]
            patch_np = np.squeeze(patch_np)
            
            if patch_np.ndim == 3 and patch_np.shape[0] in [1, 3, 4]:
                patch_np = np.transpose(patch_np, (1, 2, 0))
            if patch_np.ndim == 2:
                patch_np = np.stack((patch_np,)*3, axis=-1)
            if patch_np.shape[-1] == 4:
                patch_np = patch_np[..., :3]

            # Background pixel check
            if patch_np.mean() > BG_THRESHOLD:
                continue

            patch_pil = Image.fromarray(patch_np)
            dab_val = rgb2hed(patch_np)[..., 2].mean()
            buf.append((preprocess(patch_pil), dab_val, (left, top)))

            if len(buf) == BATCH_SIZE:
                f, d, c = flush_batch(buf)
                features_list.append(f); dab_list.extend(d); valid_coords.extend(c)
                buf = []

        if buf:
            f, d, c = flush_batch(buf)
            features_list.append(f); dab_list.extend(d); valid_coords.extend(c)

        tif.close()
        del img_array, mask  # Free up RAM immediately
        gc.collect()

        if not features_list:
            print(f"  ⚠️ No valid patches found for {slide_name} after masking.")
            continue
            
        # 3. Assemble and Cluster
        all_features = torch.cat(features_list, dim=0)
        coords_arr   = np.array(valid_coords)
        dab_arr      = np.array(dab_list)
        
        print(f"  > {all_features.shape[0]} valid patches. Running PCA & Leiden...")
        
        X = all_features.float()
        
        # PCA
        q_rank = min(128, X.shape[0] - 1)  # safe fallback if sample is tiny
        U, S, V = torch.pca_lowrank(X, q=q_rank, niter=4)
        X_pca = (X @ V).numpy()
        
        # Build AnnData
        adata = ad.AnnData(X=X.numpy())
        adata.obs["sample"] = sample_id
        adata.obs["marker"] = marker
        adata.obs["x"]      = coords_arr[:, 0].astype(float)
        adata.obs["y"]      = coords_arr[:, 1].astype(float)
        adata.obs["dab"]    = dab_arr
        adata.obsm["X_pca"] = X_pca
        
        # Scanpy workflows (Neighbors + Leiden)
        # Using a low n_neighbors fallback if patches are very few
        k_neigh = min(15, adata.n_obs - 1)
        if k_neigh > 2:
            sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=k_neigh)
            sc.tl.leiden(adata, resolution=0.5)
        else:
            adata.obs["leiden"] = "0"  # Dummy cluster if literally < 4 patches

        # Save to disk
        adata.write_h5ad(adata_out)
        print(f"  💾 Saved {adata.n_obs} patches to {adata_out}")

        # Final memory wipe per slide
        del adata, X, U, S, V, X_pca, all_features
        gc.collect()
        if torch.backends.mps.is_available():
            torch.mps.empty_cache()

print("\n✅ All markers and samples processed!")