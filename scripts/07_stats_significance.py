import glob
import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests

from config_utils import load_config, parse_config_arg, resolve_path


args = parse_config_arg("Calculate cross-sample significance for spatial enrichment matrices.")
config = load_config(args.config)

STATS_DIR = resolve_path(config, "spatial_stats_dir")
MIN_SAMPLES_FOR_PVALUE = config["significance"]["min_samples_for_pvalue"]

def calculate_proper_cohort_significance():
    print("Loading per-sample Z-score matrices...")
    csv_paths = glob.glob(os.path.join(STATS_DIR, "marker_status_zscores_sample_*.csv"))
    
    if not csv_paths:
        raise FileNotFoundError(f"No Z-score CSVs found in {STATS_DIR}.")

    # Dictionary: { ('FBXW7_high', 'MYC_low'): [2.1, 1.5, 3.0, ...], ... }
    interaction_zscores = {}
    
    for path in csv_paths:
        df = pd.read_csv(path, index_col=0)
        for row_idx in df.index:
            for col_idx in df.columns:
                str_row = str(row_idx)
                str_col = str(col_idx)

                # --- DEDUPLICATION FILTER ---
                # Because the matrix is symmetric (A to B == B to A), 
                # we strictly enforce alphabetical order. 
                # This automatically skips exact self-matches and reverse duplicates.
                if str_row >= str_col:
                    continue
                # ---------------------------------

                # --- HOMOTYPIC MARKER FILTER ---
                # Extract the base marker name (e.g., 'FBXW7' from 'FBXW7_high')
                marker_1_base = str_row.split('_')[0]
                marker_2_base = str_col.split('_')[0]
                
                # If they belong to the same base marker, skip this pair
                if marker_1_base == marker_2_base:
                    continue
                # ------------------------------------

                pair = (str_row, str_col)
                if pair not in interaction_zscores:
                    interaction_zscores[pair] = []
                
                val = df.loc[row_idx, col_idx]
                
                # Only append non-zero, non-NaN valid interactions
                if pd.notna(val) and np.isfinite(val) and val != 0.0:
                    interaction_zscores[pair].append(val)

    print("Computing cross-sample significance for unique, heterotypic interactions...")
    
    results = []
    for pair, z_array in interaction_zscores.items():
        z_array = np.array(z_array)
        k = len(z_array)
        
        if k > 0:
            mean_z = np.mean(z_array)
            std_z = np.std(z_array, ddof=1) if k > 1 else np.nan
            
            # Only calculate P-values if enough biological replicates exist
            if k >= MIN_SAMPLES_FOR_PVALUE:
                # 1-sample t-test checking if the mean Z-score is significantly != 0
                stat, p_val = ttest_1samp(z_array, popmean=0.0)
            else:
                p_val = np.nan
            
            results.append({
                "Interaction_1": pair[0],
                "Interaction_2": pair[1],
                "Mean_Enrichment": mean_z,
                "Std_Enrichment": std_z,
                "N_Samples": k,
                "P_Value": p_val
            })
            
    results_df = pd.DataFrame(results)
    
    if len(results_df) == 0:
        print("No valid unique interactions found.")
        return results_df

    # Apply FDR Correction ONLY to valid tests
    valid_p_mask = ~results_df["P_Value"].isna()
    
    if valid_p_mask.sum() > 0:
        results_df.loc[valid_p_mask, "FDR_q_value"] = multipletests(
            results_df.loc[valid_p_mask, "P_Value"], 
            alpha=0.05, 
            method="fdr_bh"
        )[1]
    else:
        results_df["FDR_q_value"] = np.nan
        print("Warning: No interactions had N >= 3. Cannot compute FDR.")
    
    # Sort by FDR, then by absolute Mean Enrichment
    results_df["Abs_Enrichment"] = results_df["Mean_Enrichment"].abs()
    results_df.sort_values(by=["FDR_q_value", "Abs_Enrichment"], ascending=[True, False], inplace=True)
    results_df.drop(columns=["Abs_Enrichment"], inplace=True)
    
    out_path = os.path.join(STATS_DIR, "Cross_Sample_Unique_Heterotypic_Significance_FDR.csv")
    results_df.to_csv(out_path, index=False)
    
    print(f"\nDone! Cleaned statistics saved to:\n{out_path}")
    return results_df

if __name__ == "__main__":
    final_stats = calculate_proper_cohort_significance()
    print("\n--- Top Unique Interactions ---")
    print(final_stats.head(15).to_string(index=False))
