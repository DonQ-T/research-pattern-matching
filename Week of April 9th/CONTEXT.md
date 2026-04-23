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

### Phase 5: Additional Similarity Metrics (DONE)
Added Fréchet, MSE, and sMAPE alongside the original three. Implementation
lives in `analysis_extended.ipynb`; results saved in `metrics_results.pkl`.
- **Fréchet** behaves nearly identically to MSE at n=4 (19-20/20 overlap)
- **MSE** is the cleanest distance metric — strong agreement with Pearson
- **sMAPE** fails on clusterThree (zero at t=3 → sMAPE_i = 2.0 for any gene)

### Phase 6: Ensemble Ranking (DONE)
Implemented in `analysis_comparison.ipynb`. Two aggregation strategies:
- Top-N vote count (how many methods place each gene in top 20)
- Average rank across methods
clusterTwo LT is the best case: 9 of top 11 genes are unanimous 6/6.
clusterFour has 0 unanimous genes — methods cannot agree on the spike.

### Phase 7: Code Review (DONE)
Verified each metric against known test vectors, audited edge cases
(constant genes, zero patterns, extreme spikes). No NaN/Inf issues found.
Notes are in `analysis_comparison.ipynb`.

### Presentation (DONE — weekly update format)
`presentation.tex` / `presentation.pdf` — restructured for weekly updates:
no methodology slides, results-focused. 3 slides per cluster (ensemble plot,
overlap heatmap, ensemble rankings table), LT data only. Appendix retains
all individual method plots. Compile with `pdflatex presentation.tex`.

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

## Recent Activity (2026-04-08)
- Pulled latest commits from GitHub (50 files: phases 5-7, math foundations,
  presentation, plots)
- Restructured `presentation.tex` to weekly-update format (removed
  methodology slides, made 3 LT slides per cluster: ensemble / heatmap /
  rankings table). Recompiled to `presentation.pdf`.
- Replied to Zihan's email with normalization explanation + GitHub repo link
- Confirmed `metrics_results.pkl` requires pandas >= 3.0 to load (use the
  `quant` conda env: `C:/Users/mjacob28/AppData/Local/miniconda3/envs/quant/python.exe`)

## Still TODO

### Phase 8: Deeper Analysis (lower priority)
- Investigate clusterFour divergence (0 overlap between any method pair) —
  needs denser temporal sampling or peak-detection approach
- Statistical significance testing (p-values for Pearson r at n=4, only 2 DOF)
- Pathway enrichment analysis on top gene lists (KEGG, GO) to validate
  biological relevance
- Cross-reference ensemble top genes with known TCR pathway genes
- Multivariate time-series similarity (Zihan is also exploring this per
  Prof. Wang's suggestion)

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
