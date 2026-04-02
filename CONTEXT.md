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

## Still TODO
- Read and incorporate findings from the time series article
- Dig deeper into differences between similarity measures
- Investigate clusterFour divergence
- Consider statistical significance testing (p-values for Pearson r at n=4)
- Possible: pathway enrichment analysis on top gene lists (KEGG, GO)

## Files
- `analysis.py` — Full script (plt.savefig, batch mode)
- `analysis.ipynb` — Jupyter notebook (inline plots, good for exploration)
- `plots/` — 32 saved PNGs (LT plots are the better ones to show)
- `SUMMARY_REPORT.md` — Detailed results summary
- `CONTEXT.md` — This file
- `CLAUDE.md` — Original implementation instructions

## People
- **Dr. Nancy Guo** — PI, flagged the NMF direction issue
- **Ethan Kyi** — Wrote the original NMF analysis (`analysis.py` original code)
