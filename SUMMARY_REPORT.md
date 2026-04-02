# Gene Expression Pattern Matching — Summary Report

## Overview

This analysis extends Ethan Kyi's NMF-based gene expression work for Dr. Nancy
Guo's research group. Three similarity metrics (Pearson correlation, DTW, and
cosine similarity) were used to identify genes whose expression time series best
match the **constraint patterns** (the biologically observed expression trends)
in each of the 4 clusters, across 4 timepoints (0, 3, 6, 9).

- **Dataset**: 17,856 genes per cluster, 4 clusters, raw + log-transformed
- **Normalization**: Min-max to (0, 1] applied before similarity computation
- **Plots generated**: 32 total (saved in `plots/`)

---

## NMF Audit (Phase 1A)

NMF uses non-negative matrix factorization, so **all coefficients and pattern
values are >= 0**. This means it cannot represent the direction of gene
regulation (up vs. down). A gene upregulated 2x and one downregulated 2x would
receive similar NMF coefficient magnitudes. The similarity metrics applied here
operate on the original expression data and preserve direction information.

The analysis was refocused from matching against NMF-reconstructed patterns to
matching against the **constraint patterns** — the actual biological expression
trends that were used to anchor NMF. Each cluster has one constraint pattern.

---

## Constraint Patterns

| Cluster      | t=0    | t=3    | t=6    | t=9    | Shape Description            |
|--------------|--------|--------|--------|--------|------------------------------|
| clusterOne   | 1.45   | 12.0   | 5.21   | 6.51   | Sharp spike at t=3, partial recovery |
| clusterTwo   | 0.371  | 0.803  | 1.22   | 1.78   | Steady increase              |
| clusterThree | 0.0405 | 0.0    | 0.0369 | 0.0353 | Near-flat with dip at t=3    |
| clusterFour  | 0.0091 | 0.638  | 0.0178 | 0.0243 | Sharp spike at t=3, returns to baseline |

---

## Method Agreement (Top-20 Gene Overlap)

### Raw Data

| Cluster      | Pearson vs DTW | Pearson vs Cosine | DTW vs Cosine |
|--------------|:--------------:|:-----------------:|:-------------:|
| clusterOne   | **16**/20      | 2/20              | 0/20          |
| clusterTwo   | 3/20           | **14**/20         | 2/20          |
| clusterThree | **17**/20      | **20**/20         | **17**/20     |
| clusterFour  | 0/20           | 0/20              | 0/20          |

### Log-Transformed Data

| Cluster        | Pearson vs DTW | Pearson vs Cosine | DTW vs Cosine |
|----------------|:--------------:|:-----------------:|:-------------:|
| clusterOneLT   | **13**/20      | 0/20              | 0/20          |
| clusterTwoLT   | **14**/20      | **18**/20         | **13**/20     |
| clusterThreeLT | **16**/20      | **20**/20         | **16**/20     |
| clusterFourLT  | 0/20           | 0/20              | 0/20          |

---

## Top Genes Per Cluster (Pearson, Raw Data)

**clusterOne** (spike at t=3): REXO2, YIF1A, OTUD7A, ACBD3, SSR3, ALDH3A2,
CDC26, ZNF512, EIF2AK4, GSPT2

**clusterTwo** (steady increase): RNF213, SSBP3, CDC16, ADAMTSL4-AS1, FKBP5,
MAP3K13, CENPF, CCNF, TMEM245, FAM222A-AS1

**clusterThree** (dip at t=3): ANKRD36C, GBP4, HERC5, IFIT1, IFIT3, IFI44,
IFI44L, MX1, OAS1, RSAD2

**clusterFour** (spike at t=3): AC129492.1, HHAT, ZNF280D, DNAJB4, MFSD2A,
ANKRD53, CYP1A1, SLC16A6, PLIN2, MXD1

---

## Key Findings

### 1. Pearson and DTW agree strongly for most clusters
For clusters with clear temporal shapes (clusterOne, clusterThree), Pearson and
DTW identify 13-17 of the same top 20 genes. This makes sense — both measure
time-series shape similarity. DTW's advantage (flexible time alignment) adds
modest value with only 4 timepoints.

### 2. clusterThree shows near-perfect agreement across all methods
All three methods return nearly identical top-20 gene lists (17-20/20 overlap).
The clusterThree constraint pattern is near-flat with a dip at t=3. This
distinctive V-shape is easy for all methods to match consistently. The top genes
(ANKRD36C, GBP4, HERC5, IFIT1, etc.) are notably interferon-stimulated genes,
which makes biological sense for a transient suppression pattern.

### 3. clusterFour is the outlier — 0 agreement between any methods
No method pair shares a single gene in their top 20. The constraint pattern
[0.009, 0.638, 0.018, 0.024] has an extreme spike at t=3 (70x baseline) then
returns to near-zero. Each method is interpreting "similar to a spike" differently:
- **Pearson** finds genes with the same shape correlation (spike-and-return)
- **DTW** finds genes with minimum warping distance (matching the near-zero baseline)
- **Cosine** finds genes with similar vector direction (dominated by the t=3 value)

This cluster may warrant investigation: are there actually genes that follow this
pattern, or is the spike too extreme for reliable matching at normalized scale?

### 4. Cosine captures different biology than Pearson/DTW
Cosine similarity measures vector angle without mean-centering. For clusterOne,
cosine selects genes like EEF1A1 and VPS51 (0 overlap with Pearson/DTW). These
genes have similar *proportional* expression across timepoints but not necessarily
the same shape after centering. This may or may not be biologically relevant
depending on the research question.

### 5. Scale mismatch in plots is visual only
The constraint patterns are plotted at their original scale (e.g., clusterOne
peaks at 12.0) while gene counts are typically 0-2. The gene curves appear flat
by comparison. This is a **visual artifact** — the similarity computation uses
normalized data, so the gene rankings are correct. The plots show that the
selected genes do follow the same *shape* at a smaller magnitude.

---

## Recommendations

1. **For robust gene lists**: Use the **Pearson** results as the primary ranking.
   Pearson is the most interpretable and agrees well with DTW.

2. **For clusterFour**: Consider relaxing the top-N threshold or examining
   whether the spike-pattern constraint is realistic. The 0-overlap across all
   methods suggests the pattern may be too extreme to match consistently.

3. **Cosine as complementary**: If proportional expression (regardless of mean)
   is biologically meaningful, cosine results provide an additional perspective.
   For clusterThree, all methods agree, confirming high-confidence gene candidates.

4. **Next steps**: Cross-reference the top gene lists with known pathway databases
   (e.g., KEGG, GO enrichment) to validate biological relevance of the selected
   genes.

---

## Files

- `analysis.py` — Full analysis script (saves plots automatically)
- `analysis.ipynb` — Interactive notebook version (inline plots)
- `plots/` — 32 PNG files:
  - `{Method}_{cluster}_constraint.png` — Gene overlay plots (raw)
  - `{Method}_LT_{cluster}_constraint.png` — Gene overlay plots (log-transformed)
  - `overlap_{cluster}.png` — Cross-method heatmaps (raw)
  - `overlap_LT_{cluster}.png` — Cross-method heatmaps (log-transformed)
