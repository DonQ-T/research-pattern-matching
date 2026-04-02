# Research: Gene Expression Pattern Matching

## Context
Extending Ethan Kyi's NMF-based gene expression analysis with additional
similarity metrics. This is for Dr. Nancy Guo's research group at Binghamton.
The goal is to evaluate whether Pearson correlation, DTW, and cosine similarity
capture biologically meaningful patterns better than (or complementary to) NMF.

## Data Structure
- 4 clusters of ~17,857 genes across 4 timepoints (0, 3, 6, 9)
- `analysis data/gene_counts/` — raw expression counts per cluster (4 CSV files)
- `analysis data/gene_countsLT/` — log-transformed counts (4 CSV files)
- `analysis data/NMF run results/clusters/` — NMF output (patterns + per-gene coefficients)
- `analysis data/NMF run results/clustersLT/` — NMF output for log-transformed data
- `analysis data/gene_list.csv` — full gene list (~17,857 genes)
- `analysis.py` — Ethan's existing analysis code

CSV format: rows are genes (gene name as index), columns are 4 timepoints.
NMF results: header lines starting with `#` contain pattern vectors, `row-N` lines contain per-gene coefficients.

## Required Modifications (in order)

### Phase 1: Audit and prep

#### 1A. Audit NMF fold-change handling
- Check whether the current NMF approach uses absolute fold change
  (losing up/down-regulation direction as Nancy flagged)
- Document findings as comments in the code
- If direction is lost, note which genes/patterns are affected

#### 1B. Fix file paths
- Ethan's paths are relative to a different directory structure
- Update all paths to work from this project root (`analysis data/...`)

#### 1C. Add normalization to (0, 1]
- Implement min-max normalization for gene expression time series
- Apply to both raw and LT data before similarity comparison
- Preserve unnormalized data for reference

**Validation:** After Phase 1, print shape and summary stats of loaded data
for one cluster. Confirm gene counts match expected ~17,857. Confirm
normalized values fall in (0, 1]. Run this check before proceeding.

### Phase 2: Implement similarity metrics

#### 2A. Pearson correlation
- For each pattern in each cluster, compute Pearson correlation between
  the pattern time series and every gene's expression time series
- Rank genes by correlation (highest = most similar)
- Store results in same structure as Ethan's top_genes dicts

#### 2B. DTW (Dynamic Time Warping)
- Use `dtw-python` or `tslearn` library
- Normalize both pattern and gene series before computing DTW
- Compute DTW distance between each pattern and each gene's time series
- Rank genes by DTW distance (lower = more similar)

#### 2C. Cosine similarity
- Compute cosine similarity between pattern vectors and gene expression vectors
- Rank genes (highest similarity = most similar)

**Validation:** After Phase 2, for clusterOne/pattern1, print the top 10 genes
from each method (NMF, Pearson, DTW, Cosine) side by side. Sanity check:
there should be *some* overlap between methods — if a method returns
completely disjoint results, investigate before continuing.

### Phase 3: Visualization and comparison

#### 3A. Per-method plots
- Generate plots matching Ethan's style (matplotlib, same color scheme)
- For each metric, plot top-N gene expression curves overlaid with the pattern
- Use x = [0, 3, 6, 9] timepoints consistently

#### 3B. Cross-method comparison
- For each cluster/pattern, produce a summary showing how top-N gene
  rankings differ across methods (NMF vs Correlation vs DTW vs Cosine)
- Overlap counts (e.g., "15 of top 20 shared between Pearson and NMF")
- Consider a heatmap or Venn-style visualization of ranking overlap

**Validation:** After Phase 3, visually confirm at least one plot per method
renders correctly. Check axis labels, legends, and that gene curves
actually overlay the pattern (not a blank or misaligned plot).

### Phase 4: Final validation
- Run the complete script end-to-end
- Confirm all plots generate without errors
- Print a final summary table: for each cluster and pattern, show the
  number of overlapping genes in the top 20 across all method pairs

## Execution Strategy — Use Subagents

When implementing this plan, use subagents to parallelize independent work:

- **Phase 2 parallelism:** Implement Pearson (2A), DTW (2B), and Cosine (2C)
  as independent subagents since they share no dependencies beyond the
  normalized data from Phase 1. Each subagent should write its own function(s)
  and return the results structure.
- **Phase 3 parallelism:** Plot generation for different methods/clusters
  can be parallelized across subagents.
- **Research subagents:** If you need to look up DTW library APIs or
  normalization best practices, spawn a research subagent rather than
  blocking the main flow.

Do NOT use subagents for sequential/dependent work (e.g., Phase 1 must
complete before Phase 2 begins).

## Constraints
- **Preserve Ethan's existing code** — add new functions below his, don't rewrite
- Keep the same plotting style (matplotlib) for visual consistency
- Use the same timepoint convention: x = [0, 3, 6, 9]
- All new functions should work on both raw and LT data (pass `lt=False` flag)
- Print progress messages so the user can track execution
- Install any needed packages (dtw-python, tslearn, scipy) at the start
