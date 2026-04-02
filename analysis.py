import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- Phase 1A: NMF fold-change audit ---
# NMF (Non-negative Matrix Factorization) decomposes the expression matrix V
# into two non-negative matrices W (coefficients) and H (patterns): V ≈ W * H.
# Because both W and H are constrained to be non-negative, NMF CANNOT represent
# the direction of regulation (up vs. down). A gene that is upregulated 2x and
# one that is downregulated 2x would receive similar coefficient magnitudes.
# This effectively uses absolute fold change, losing sign information.
# Nancy flagged this — the new similarity metrics (Pearson, DTW, Cosine) applied
# to the original expression data will preserve direction and may capture
# biologically meaningful patterns that NMF misses.

# --- Phase 1B: Fix file paths (relative to project root) ---
results = "analysis data/NMF run results/clusters"
resultsLT = "analysis data/NMF run results/clustersLT"
genes_csv = "analysis data/gene_list.csv"

genes = pd.read_csv(genes_csv, header=None).squeeze().tolist()

def extract_coeff(folder) -> pd.DataFrame:
    cluster_dfs = {}
    
    for filename in sorted(os.listdir(folder)):
        path = os.path.join(folder, filename)
        rows = []
        
        with open(path) as file:
            for line in file:
                if line.startswith("row-"):
                    values = [float(x) for x in line[4:].strip().split(",")[1:4]]
                    rows.append(values)
        
        num_patterns = len(rows[0]) if rows else 0
        pattern_cols = [f"pattern{i+1}" for i in range(num_patterns)]
        
        cluster_name = os.path.splitext(filename)[0].split('_')[0]
        cluster_dfs[cluster_name] = pd.DataFrame(rows, columns=pattern_cols, index=genes)
        cluster_dfs[cluster_name].index.name = "gene"
    
    return cluster_dfs

cluster_dfs = extract_coeff(results)
clusterLT_dfs = extract_coeff(resultsLT)
#print(cluster_dfs)

def extract_patterns(folder):
    patterns_df = {}
    
    for filename in sorted(os.listdir(folder)):
        path = os.path.join(folder, filename)
        cluster_name = os.path.splitext(filename)[0].split('_')[0]
        
        if cluster_name not in patterns_df:
            patterns_df[cluster_name] = {}
        
        with open(path) as file:
            for line in file:
                if line.startswith('#"'):
                    parts = line.strip().split(',')
                    parts = parts[1].split()
                    vals = [float(x) for x in parts]
                    pattern_key = f"pattern{len(patterns_df[cluster_name]) + 1}"
                    patterns_df[cluster_name][pattern_key] = vals
    
    return patterns_df

patterns = extract_patterns(results)
patternsLT = extract_patterns(resultsLT)
#print(patterns)
#print(patterns)

def top_genes_per_pattern(cluster_df, n=20):
    top_ranked_genes = {}

    for cluster_name, df in cluster_df.items():
        top_ranked_genes[cluster_name] = {}
        for pattern in df.columns:
            top_ranked_genes[cluster_name][pattern] = (df[pattern].sort_values(ascending = False).head(n))
    return top_ranked_genes


top_genes = top_genes_per_pattern(cluster_dfs)
top_genesLT = top_genes_per_pattern(clusterLT_dfs)

#print(top_genesLT)
#print(top_genes["clusterOne"])

#plot patterns identified in each cluster

constraints = {
    "clusterOne": [1.45,12.0,5.21,6.51],
    "clusterTwo": [0.371,0.803,1.22,1.78],
    "clusterThree": [0.0405,0,0.0369,0.0353],
    "clusterFour": [0.00913,0.638,0.0178,0.0243]
}

constraintsLT = {
    "clusterOneLT": [.894475,2.568650,1.826,2.016452],
    "clusterTwoLT": [0.036407,0.589570,0.798347,1.023510],
    "clusterThreeLT": [0.039660,0.0,0.036266,0.034714],
    "clusterFourLT": [0.009090,0.493705,0.017601,0.024007]
}

def plot_patterns(pattern_df, lt=False):
    x = [0, 3, 6, 9]
    constraint = {}
    constraint = constraintsLT if lt else constraints
    for cluster in pattern_df:
        print(pattern_df[cluster])
        for pattern in pattern_df[cluster]:
            print(pattern_df[cluster][pattern])
            plt.plot(x, pattern_df[cluster][pattern], label = f"{pattern}")
        print('\n')
        plt.plot(x, constraint[cluster], color="crimson", linewidth=2.5, linestyle="--", label=f"{cluster} constraint")
            
        plt.title(f"{cluster} patterns")
        plt.xlabel("timepoint")
        plt.ylabel("coefficients")
        plt.legend()
        plt.xticks(x)
        plt.show()

#plot_patterns(patterns)
#plot_patterns(patternsLT, True)


def plot_genes_per_pattern(patterns, top_genes, lt=False, n=20):
    x = [0, 3, 6, 9]
    suffix = "LT" if lt else ""
    gene_folder = f"analysis data/gene_counts{suffix}/"

    for cluster_name, pattern_dict in patterns.items():
        print(cluster_name)
        gene_file = os.path.join(gene_folder, f"{cluster_name}_annotated.csv")
        gene_counts = pd.read_csv(gene_file, index_col=0)
        gene_counts.columns = x
        

        for pattern_name, pattern_vals in pattern_dict.items():
            top = top_genes[cluster_name][pattern_name].head(n)
            print(top)
            top_gene_names = top.index.tolist()

            fig, ax = plt.subplots(figsize=(10, 6))

            for gene in top_gene_names:
                if gene in gene_counts.index:
                    ax.plot(x, gene_counts.loc[gene], color="steelblue",
                            alpha=0.4, linewidth=1, label=gene)

            #pattern_vals = [5 * x for x in pattern_vals]
            ax.plot(x, pattern_vals, color="crimson", linewidth=2.5, linestyle="--", label=f"{pattern_name}")
            
            ax.set_title(f"{cluster_name} — {pattern_name} | Top {n} genes {'(LT)' if lt else ''}")
            ax.set_xlabel("Timepoint")
            ax.set_ylabel("Gene counts")
            ax.set_xticks(x)

            handles, labels = ax.get_legend_handles_labels()
            pattern_handle = [(h, l) for h, l in zip(handles, labels) if "scaled" in l]
            gene_handles = [(h, l) for h, l in zip(handles, labels) if "scaled" not in l]

            ax.legend(
                [ph[0] for ph in pattern_handle] + [gene_handles[0][0]],
                [ph[1] for ph in pattern_handle] + [f"Top {n} genes (n={len(gene_handles)})"],
                loc="upper right"
            )

            plt.tight_layout()
            plt.show()

#plot_genes_per_pattern(patterns, top_genes)
#plot_genes_per_pattern(patternsLT, top_genesLT, lt=True)


# ============================================================================
# Phase 1C: Load gene expression data and normalize to (0, 1]
# ============================================================================

def load_gene_counts(lt=False):
    """Load gene expression counts for all clusters. Returns dict of DataFrames."""
    suffix = "LT" if lt else ""
    folder = f"analysis data/gene_counts{suffix}/"
    x = [0, 3, 6, 9]
    cluster_counts = {}
    for filename in sorted(os.listdir(folder)):
        if not filename.endswith(".csv"):
            continue
        path = os.path.join(folder, filename)
        df = pd.read_csv(path, index_col=0)
        df.columns = x
        cluster_name = filename.split("_")[0]
        cluster_counts[cluster_name] = df
    return cluster_counts


def normalize_01(cluster_counts):
    """Min-max normalize each gene's time series to (0, 1].

    For each gene row: normalized = (x - min + eps) / (max - min + eps)
    This ensures all values are in (0, 1] — strictly positive, max maps to 1.
    Genes with constant expression across all timepoints get value 1.0.
    """
    eps = 1e-10
    normalized = {}
    for cluster_name, df in cluster_counts.items():
        row_min = df.min(axis=1)
        row_max = df.max(axis=1)
        row_range = row_max - row_min
        # For constant rows (range=0), set to 1.0 via eps cancellation
        norm_df = df.sub(row_min, axis=0).add(eps).div(row_range.add(eps), axis=0)
        # Genes where all values are identical → set to 1.0
        constant_mask = row_range == 0
        norm_df.loc[constant_mask] = 1.0
        normalized[cluster_name] = norm_df
    return normalized


# Load raw and LT expression data
print("=" * 60)
print("Phase 1: Loading and normalizing gene expression data")
print("=" * 60)

gene_counts = load_gene_counts(lt=False)
gene_countsLT = load_gene_counts(lt=True)

# Normalize
gene_counts_norm = normalize_01(gene_counts)
gene_countsLT_norm = normalize_01(gene_countsLT)

# --- Phase 1 Validation ---
print("\n--- Phase 1 Validation ---")
sample_cluster = list(gene_counts.keys())[0]
sample_df = gene_counts[sample_cluster]
sample_norm = gene_counts_norm[sample_cluster]
print(f"Cluster: {sample_cluster}")
print(f"  Shape: {sample_df.shape} (expected ~17857 genes x 4 timepoints)")
print(f"  Gene count: {sample_df.shape[0]}")
print(f"  Raw data summary:\n{sample_df.describe()}\n")
print(f"  Normalized data summary:\n{sample_norm.describe()}")
print(f"  Normalized min: {sample_norm.min().min():.10f} (should be > 0)")
print(f"  Normalized max: {sample_norm.max().max():.10f} (should be <= 1)")
assert sample_norm.min().min() > 0, "Normalized values must be > 0"
assert sample_norm.max().max() <= 1.0, "Normalized values must be <= 1"
print("  [OK] Normalization validated: all values in (0, 1]")
print(f"  Total genes across all clusters: {sum(df.shape[0] for df in gene_counts.values())}")
print()


# ============================================================================
# Phase 2: Similarity metrics
# ============================================================================

from scipy import stats as scipy_stats
from scipy.spatial.distance import cosine as cosine_distance
from dtw import dtw as dtw_func


# --- Phase 2A: Pearson Correlation ---
# Computes Pearson r between each NMF pattern and every gene's normalized
# expression. Preserves up/down-regulation direction (unlike NMF).

def top_genes_pearson(patterns, gene_counts_norm, n=20):
    results = {}
    for cluster_name, cluster_patterns in patterns.items():
        print(f"  [Pearson] Processing {cluster_name}")
        df = gene_counts_norm[cluster_name]
        gene_matrix = df.values
        gene_names = df.index.tolist()
        n_genes = gene_matrix.shape[0]
        cluster_results = {}

        for pattern_name, pattern_vec in cluster_patterns.items():
            pattern_arr = np.array(pattern_vec, dtype=float)
            if np.std(pattern_arr) == 0:
                r_values = np.zeros(n_genes)
            else:
                r_values = np.empty(n_genes)
                for i in range(n_genes):
                    gene_vec = gene_matrix[i]
                    if np.std(gene_vec) == 0:
                        r_values[i] = 0.0
                    else:
                        r, _ = scipy_stats.pearsonr(pattern_arr, gene_vec)
                        r_values[i] = r if np.isfinite(r) else 0.0

            r_series = pd.Series(r_values, index=gene_names, name=pattern_name)
            top_n = r_series.sort_values(ascending=False).head(n)
            cluster_results[pattern_name] = top_n
            print(f"    {pattern_name}: top = {top_n.index[0]} (r={top_n.iloc[0]:.4f})")

        results[cluster_name] = cluster_results
    return results


# --- Phase 2B: DTW (Dynamic Time Warping) ---
# DTW distance between normalized pattern and gene time series.
# Lower distance = more similar.

def top_genes_dtw(patterns, gene_counts_norm, n=20):
    results = {}
    eps = 1e-10
    for cluster_name, cluster_patterns in patterns.items():
        print(f"  [DTW] Processing {cluster_name}")
        results[cluster_name] = {}
        df = gene_counts_norm[cluster_name]
        gene_matrix = df.values

        for pattern_name, pattern_vals in cluster_patterns.items():
            p = np.array(pattern_vals, dtype=float)
            p_norm = (p - p.min() + eps) / (p.max() - p.min() + eps)

            distances = np.empty(len(df), dtype=float)
            for i in range(len(df)):
                gene_vec = gene_matrix[i].astype(float)
                alignment = dtw_func(p_norm, gene_vec)
                distances[i] = alignment.distance

            dist_series = pd.Series(distances, index=df.index)
            top_n = dist_series.nsmallest(n)
            results[cluster_name][pattern_name] = top_n
            print(f"    {pattern_name}: top = {top_n.index[0]} (dist={top_n.iloc[0]:.4f})")

    return results


# --- Phase 2C: Cosine Similarity ---
# Cosine similarity between pattern vectors and gene expression vectors.
# Higher similarity = more similar.

def top_genes_cosine(patterns, gene_counts_norm, n=20):
    results = {}
    for cluster_name, cluster_patterns in patterns.items():
        print(f"  [Cosine] Processing {cluster_name}")
        results[cluster_name] = {}
        df = gene_counts_norm[cluster_name]
        gene_matrix = df.values
        gene_names = df.index.tolist()

        for pattern_name, pattern_vec in cluster_patterns.items():
            pattern_arr = np.array(pattern_vec, dtype=float)
            if np.linalg.norm(pattern_arr) == 0:
                similarities = np.zeros(len(gene_names))
            else:
                similarities = np.empty(len(gene_names), dtype=float)
                for i, gene_vec in enumerate(gene_matrix):
                    if np.linalg.norm(gene_vec) == 0:
                        similarities[i] = 0.0
                    else:
                        similarities[i] = 1.0 - cosine_distance(pattern_arr, gene_vec)

            sim_series = pd.Series(similarities, index=gene_names)
            top_n = sim_series.nlargest(n)
            results[cluster_name][pattern_name] = top_n
            print(f"    {pattern_name}: top = {top_n.index[0]} (sim={top_n.iloc[0]:.4f})")

    return results


# --- Convert constraints to pattern dict format ---
# The constraints ARE the biological patterns we want to match genes against.
# One pattern per cluster (the observed expression trend).
constraint_patterns = {k: {"constraint": v} for k, v in constraints.items()}
constraint_patternsLT = {k: {"constraint": v} for k, v in constraintsLT.items()}

# --- Run Phase 2 ---
print("=" * 60)
print("Phase 2: Computing similarity metrics (against constraint patterns)")
print("=" * 60)

print("\n[Pearson correlation]")
pearson_top = top_genes_pearson(constraint_patterns, gene_counts_norm, n=20)
pearson_topLT = top_genes_pearson(constraint_patternsLT, gene_countsLT_norm, n=20)

print("\n[DTW distance]")
dtw_top = top_genes_dtw(constraint_patterns, gene_counts_norm, n=20)
dtw_topLT = top_genes_dtw(constraint_patternsLT, gene_countsLT_norm, n=20)

print("\n[Cosine similarity]")
cosine_top = top_genes_cosine(constraint_patterns, gene_counts_norm, n=20)
cosine_topLT = top_genes_cosine(constraint_patternsLT, gene_countsLT_norm, n=20)

# --- Phase 2 Validation ---
print("\n--- Phase 2 Validation ---")
val_cluster = list(constraint_patterns.keys())[0]
print(f"Comparing top 10 genes for {val_cluster}/constraint:\n")

pears_top10 = pearson_top[val_cluster]["constraint"].head(10).index.tolist()
dtw_top10 = dtw_top[val_cluster]["constraint"].head(10).index.tolist()
cos_top10 = cosine_top[val_cluster]["constraint"].head(10).index.tolist()

print(f"{'Rank':<5} {'Pearson':<20} {'DTW':<20} {'Cosine':<20}")
print("-" * 65)
for i in range(10):
    print(f"{i+1:<5} {pears_top10[i]:<20} {dtw_top10[i]:<20} {cos_top10[i]:<20}")

# Check overlap between methods
all_methods = {"Pearson": set(pears_top10), "DTW": set(dtw_top10),
               "Cosine": set(cos_top10)}
print(f"\nOverlap (top 10):")
method_names = list(all_methods.keys())
for i in range(len(method_names)):
    for j in range(i + 1, len(method_names)):
        m1, m2 = method_names[i], method_names[j]
        overlap = len(all_methods[m1] & all_methods[m2])
        print(f"  {m1} vs {m2}: {overlap}/10 shared")
print()


# ============================================================================
# Phase 3: Visualization and cross-method comparison
# ============================================================================

def plot_genes_per_method(method_name, method_top, patterns, gene_counts,
                          lt=False, n=20):
    """Plot top-N gene expression curves overlaid with the pattern for a given method.

    Matches Ethan's plotting style (matplotlib, steelblue genes, crimson pattern).
    Saves figures to plots/ directory instead of plt.show() for batch execution.
    """
    x = [0, 3, 6, 9]
    suffix = "LT" if lt else ""
    tag = f"{method_name}{'_LT' if lt else ''}"
    os.makedirs("plots", exist_ok=True)

    for cluster_name, pattern_dict in patterns.items():
        gene_file = f"analysis data/gene_counts{suffix}/{cluster_name}_annotated.csv"
        gene_df = pd.read_csv(gene_file, index_col=0)
        gene_df.columns = x

        for pattern_name, pattern_vals in pattern_dict.items():
            top_series = method_top[cluster_name][pattern_name].head(n)
            top_gene_names = top_series.index.tolist()

            fig, ax = plt.subplots(figsize=(10, 6))

            plotted = 0
            for gene in top_gene_names:
                if gene in gene_df.index:
                    ax.plot(x, gene_df.loc[gene], color="steelblue",
                            alpha=0.4, linewidth=1)
                    plotted += 1

            ax.plot(x, pattern_vals, color="crimson", linewidth=2.5,
                    linestyle="--", label=f"{pattern_name} (pattern)")

            ax.set_title(f"{method_name} | {cluster_name} - {pattern_name} | "
                         f"Top {n} genes {'(LT)' if lt else ''}")
            ax.set_xlabel("Timepoint")
            ax.set_ylabel("Gene counts")
            ax.set_xticks(x)

            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color="crimson", linewidth=2.5, linestyle="--",
                       label=f"{pattern_name}"),
                Line2D([0], [0], color="steelblue", alpha=0.4,
                       label=f"Top {n} genes (n={plotted})")
            ]
            ax.legend(handles=legend_elements, loc="upper right")

            plt.tight_layout()
            fname = f"plots/{tag}_{cluster_name}_{pattern_name}.png"
            plt.savefig(fname, dpi=150)
            plt.close(fig)
            print(f"  Saved: {fname}")


def plot_cross_method_heatmap(patterns, pearson_top, dtw_top,
                              cosine_top, lt=False, n=20):
    """For each cluster, create a heatmap of pairwise overlap counts
    across methods (Pearson, DTW, Cosine) for the top-N genes."""
    method_names = ["Pearson", "DTW", "Cosine"]
    method_dicts = [pearson_top, dtw_top, cosine_top]
    nm = len(method_names)

    for cluster_name, pattern_dict in patterns.items():
        for pattern_name in pattern_dict:
            gene_sets = []
            for md in method_dicts:
                genes = set(md[cluster_name][pattern_name].head(n).index.tolist())
                gene_sets.append(genes)

            overlap_matrix = np.zeros((nm, nm), dtype=int)
            for i in range(nm):
                for j in range(nm):
                    overlap_matrix[i, j] = len(gene_sets[i] & gene_sets[j])

            fig, ax = plt.subplots(figsize=(5, 4))
            im = ax.imshow(overlap_matrix, cmap="YlOrRd", vmin=0, vmax=n)
            ax.set_xticks(range(nm))
            ax.set_yticks(range(nm))
            ax.set_xticklabels(method_names)
            ax.set_yticklabels(method_names)
            ax.set_title(f"Top-{n} Gene Overlap | {cluster_name}"
                         f"{' (LT)' if lt else ''}")

            for i in range(nm):
                for j in range(nm):
                    ax.text(j, i, str(overlap_matrix[i, j]),
                            ha="center", va="center", fontsize=14,
                            color="white" if overlap_matrix[i, j] > n / 2 else "black")

            plt.colorbar(im, ax=ax, label="Shared genes")
            plt.tight_layout()
            suffix_str = "_LT" if lt else ""
            fname = f"plots/overlap{suffix_str}_{cluster_name}.png"
            plt.savefig(fname, dpi=150)
            plt.close(fig)
            print(f"  Saved: {fname}")


# --- Run Phase 3 ---
print("=" * 60)
print("Phase 3: Generating plots")
print("=" * 60)

print("\n[Per-method plots - Raw data]")
plot_genes_per_method("Pearson", pearson_top, constraint_patterns, gene_counts, lt=False)
plot_genes_per_method("DTW", dtw_top, constraint_patterns, gene_counts, lt=False)
plot_genes_per_method("Cosine", cosine_top, constraint_patterns, gene_counts, lt=False)

print("\n[Per-method plots - Log-transformed data]")
plot_genes_per_method("Pearson", pearson_topLT, constraint_patternsLT, gene_countsLT, lt=True)
plot_genes_per_method("DTW", dtw_topLT, constraint_patternsLT, gene_countsLT, lt=True)
plot_genes_per_method("Cosine", cosine_topLT, constraint_patternsLT, gene_countsLT, lt=True)

print("\n[Cross-method overlap heatmaps]")
plot_cross_method_heatmap(constraint_patterns, pearson_top, dtw_top, cosine_top,
                          lt=False)
plot_cross_method_heatmap(constraint_patternsLT, pearson_topLT, dtw_topLT,
                          cosine_topLT, lt=True)

print("\n[OK] Phase 3 complete.")


# ============================================================================
# Phase 4: Final validation — summary table
# ============================================================================

print("\n" + "=" * 60)
print("Phase 4: Final summary — top-20 overlap across all method pairs")
print("=" * 60)

method_pairs = [("Pearson", "DTW"), ("Pearson", "Cosine"), ("DTW", "Cosine")]

for label, pat_dict, pt, dt, ct in [
    ("Raw", constraint_patterns, pearson_top, dtw_top, cosine_top),
    ("LT", constraint_patternsLT, pearson_topLT, dtw_topLT, cosine_topLT)
]:
    print(f"\n--- {label} data ---")
    method_map = {"Pearson": pt, "DTW": dt, "Cosine": ct}

    pair_headers = [f"{a[:3]}v{b[:3]}" for a, b in method_pairs]
    header = f"{'Cluster':<20} " + " ".join(f"{h:>10}" for h in pair_headers)
    print(header)
    print("-" * len(header))

    for cluster_name, pattern_dict in pat_dict.items():
        overlaps = []
        for m1_name, m2_name in method_pairs:
            s1 = set(method_map[m1_name][cluster_name]["constraint"].head(20).index)
            s2 = set(method_map[m2_name][cluster_name]["constraint"].head(20).index)
            overlaps.append(len(s1 & s2))
        vals = " ".join(f"{o:>10}" for o in overlaps)
        print(f"{cluster_name:<20} {vals}")

print("\n[OK] All phases complete.")
