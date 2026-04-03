# Project Context — Gene Expression Pattern Matching

## What This Project Is
Extending Ethan Kyi's NMF-based gene expression analysis for Dr. Nancy Guo's
research group at Binghamton. The goal is to evaluate whether Pearson correlation,
DTW, and cosine similarity can identify genes whose expression matches observed
biological patterns (constraints) better than NMF alone.

## The Data
- **Source**: TCR (T-cell receptor) gene expression data
- **Structure**: ~17,856 genes across 4 timepoints (0, 3, 6, 9)
- **4 clusters**, each with raw and log-transformed (LT) versions
  - `analysis data/gene_counts/` — raw expression CSVs
  - `analysis data/gene_countsLT/` — log-transformed CSVs
  - `analysis data/NMF run results/` — Ethan's NMF output
  - `analysis data/gene_list.csv` — full gene list
- **Constraint patterns** (one per cluster) — the biologically observed
  expression trends used to anchor NMF. These are what we match genes against.

## What Was Done

### Phase 1: Audit and prep
- Confirmed NMF loses direction of regulation (up/down) because all factors
  are non-negative. It ranks genes by coefficient magnitude, not shape match.
- Fixed Ethan's file paths to work from project root
- Applied min-max normalization to (0, 1] on gene expression time series

### Phase 2: Similarity metrics
Three methods compare each gene's normalized expression to the constraint pattern:
- **Pearson correlation** — mean-centers both vectors, measures shape similarity.
  Best method for this data. Directly captures "does the gene go up/down at the
  same timepoints as the pattern?"
- **DTW (Dynamic Time Warping)** — finds optimal time alignment, measures distance.
  Largely redundant with Pearson at only 4 timepoints (13-17/20 overlap).
  Would add more value with denser time series (15+ points).
- **Cosine similarity** — measures vector angle WITHOUT mean-centering. Conflates
  high baseline expression with shape similarity. Weakest method for this use case.
  Finds genes like EEF1A1 (housekeeping) that don't actually match the shape.

### Phase 3: Visualization
- Per-method plots: top 20 genes overlaid on constraint pattern
- Cross-method overlap heatmaps
- LT plots look much better than raw (log transform compresses the scale
  so gene curves and pattern are visually comparable)

### Phase 4: Summary
- clusterThree: all methods agree (17-20/20 overlap), interferon-stimulated genes
- clusterOne/Two: Pearson-DTW agree well (13-16/20)
- clusterFour: 0 overlap across all methods — extreme spike pattern is hard to match

## Key Insight: Why NMF Was Insufficient
NMF uses absolute fold change — a gene upregulated 2x and downregulated 2x get
the same coefficient. Its top genes (RPS18, MALAT1, RPL13) are highly expressed
housekeeping genes, not genes matching the pattern shape. The new similarity
metrics compare directly against the constraint patterns and preserve direction.

## Known Limitations
- **4 timepoints is very few** — high correlations can occur by chance (only 2
  degrees of freedom for Pearson). Statistical significance of individual gene
  rankings should be interpreted carefully.
- **DTW adds little at this resolution** — essentially reduces to point-by-point
  comparison with 4 points. More valuable with 15+ timepoints.
- **Cosine doesn't mean-center** — biased toward high-expression genes regardless
  of temporal dynamics.
- **Scale mismatch in raw plots** — constraint pattern plotted at original scale
  while genes are at expression count scale. LT plots are more visually informative.

## People
- **Dr. Nancy Guo** — PI, AI/ML SUNY Empire Innovation Professor, School of Computing, Binghamton
- **Minjie Wang** — Co-advisor/collaborator, suggested similarity metrics and ensemble methods
- **Ethan Kyi** — Wrote original NMF analysis code and provided data via Google Drive/GitHub
- **Zihan** — Team member (on email thread)
- **Matt Jacob** — Analysis developer, implementing the new similarity metrics

## Still TODO

### Phase 5: Additional Similarity Metrics

#### 5A. Fréchet Distance
- Minjie recommended this as a 4th distance metric (email Mar 26)
- Measures similarity between curves considering both position and ordering
- Unlike DTW, Fréchet distance preserves temporal ordering (no warping)
- Library: `similaritymeasures` or `frechetdist` in Python
- Note: with only 4 timepoints, may behave similarly to Euclidean distance

#### 5B. MSE (Mean Squared Error)
- Nancy and Minjie discussed as a similarity measure (email Mar 28)
- Point-by-point squared difference between normalized pattern and gene series
- Simple and interpretable; lower = more similar
- Already available in sklearn (`mean_squared_error`) or trivial to implement

#### 5C. MAPE (Mean Absolute Percentage Error)
- Minjie recommended (email Mar 29)
- Measures percentage deviation at each timepoint
- Caution: undefined when true value = 0 (clusterThree has a 0 at t=3)
- Need to handle zero-value timepoints (skip or use symmetric MAPE)

### Phase 6: Ensemble Ranking Aggregation
- Nancy asked: "can we develop ensemble methods to get the best results from each method?"
- Minjie proposed two aggregation strategies:
  1. **Top-N vote count**: Count how many methods place each gene in the top N
  2. **Average/median rank**: Compute the mean or median rank of each gene across methods
- Additional ideas from Minjie: investigate statistical ensemble methods
- Output: a consensus "ensemble top-N" gene list per cluster/pattern
- This is likely the most impactful deliverable — synthesizes all individual methods

### Phase 7: Deeper Analysis
- Investigate clusterFour divergence (0 overlap between any method pair)
- Statistical significance testing (p-values for Pearson r at n=4, only 2 DOF)
- Pathway enrichment analysis on top gene lists (KEGG, GO) to validate biological relevance
- Cross-reference ensemble top genes with known TCR pathway genes

### Research References
- **SimilarityTS paper**: https://www.sciencedirect.com/science/article/pii/S2352711023002236
  - Python toolkit for standardized evaluation of time series similarity
  - Section 3.2.1 covers univariate metrics (behind paywall — access via Binghamton library)
  - Covers: KL divergence, DTW, PCA/t-SNE visualization comparisons
- **Medium post (L. Bader)**: https://medium.com/@lxbader/how-can-we-quantify-similarity-between-time-series-ed1d0b633ca0
  - Practical comparison of Euclidean, Pearson, DTW, and Compression-Based Dissimilarity
  - Key insight: Pearson best for shape matching (ignores amplitude), DTW for time-shifted
  - With 4 timepoints, DTW warping flexibility is very limited
  - Also covers SAX-based compression similarity (lower priority, more exotic)

## Files
- `analysis.py` — Full script (plt.savefig, batch mode)
- `analysis.ipynb` — Jupyter notebook (inline plots, good for exploration)
- `plots/` — 32 saved PNGs (LT plots are the better ones to show)
- `SUMMARY_REPORT.md` — Detailed results summary
- `CONTEXT.md` — This file
- `CLAUDE.md` — Implementation instructions
