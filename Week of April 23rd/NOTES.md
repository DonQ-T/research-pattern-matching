# Week of April 23 — Study Notes

Study guide for the shape vs. magnitude framing. Scan the TL;DR, then skim the rest; drill into Q&A at the end.

---

## TL;DR

1. April 16 finding: **raw MAE** (unnormalized LT distance) gave the best *visual* fit for clusters Two, Three, and Four. The 6 prior methods all run on (0,1]-normalized data, which erases magnitude — so they find the right shape at arbitrary expression levels. ClusterOne is the exception.
2. That forces a biological question: **is "similar" about direction of change (shape) or absolute expression level (magnitude)?**
3. This week's deck doesn't answer it — it lays out the two paths (Path A: shape-based gated ensemble / Path B: raw MAE) per cluster and asks the group.
4. At only 4 timepoints, "shape" = 3 transitions with 2³ = 8 possible sign classes. That's a weak filter for ~17,856 genes.

---

## Core concepts

### Normalization erases magnitude
All 6 original methods (Pearson, DTW, Cosine, Fréchet, MSE, sMAPE) are applied *after* each gene and each pattern is min-max normalized to (0,1]. That step strips the absolute expression level before any comparison happens. After normalization, two genes with baselines differing by 100× look identical — only the relative up-and-down shape survives.

### Shape classes (sign patterns)
With 4 timepoints, each gene has 3 transitions: t0→t3, t3→t6, t6→t9. Collapse each to its sign (+/−/0) and there are at most 3³ = 27 "classes" — though only a handful are populated (genes with mixed 0-transitions are rare). A shape-based metric sorts all genes into these classes with fine-grained tiebreaking.

**Percent of genes sharing the constraint's sign class**:
| Cluster | Sign class | Matching / Total | % |
|---|---|---|---|
| clusterOne | (+,−,+) | 1,530 / 17,856 | 8.6% |
| clusterTwo | (+,+,+) | 160 / 17,856 | 0.9% |
| clusterThree | (−,+,−) | 6,326 / 17,856 | **35.4%** |
| clusterFour | (+,−,+) | 1,315 / 17,856 | 7.4% |

The **35.4% on clusterThree** is the damning number — any shape-based method has ~6,300 genes that trivially match the sign class before tiebreaking. It's a weak filter.

### Which methods are inherently scale-invariant?
- **Pearson**: yes (subtracts mean, divides by std).
- **Cosine**: yes (angle-based).
- **DTW, Fréchet, MSE, sMAPE**: no — these are distance metrics. They only behaved as shape-based in the April 9 pipeline because normalization happened upstream. The appendix shows what they actually do on raw data.

### The Gate C-norm cutoff (22.946)
**What it measures.** For each method × cluster pair, "relative centroid error":
1. Compute the mean of the method's top-20 gene expression curves → a centroid curve.
2. MSE(centroid, constraint pattern) → how far the centroid sits from the pattern.
3. Divide by (pattern range)² so clusters with different expression scales are comparable.

**Formula**: `Relative Centroid Error = MSE(centroid_top20, pattern) / (max(pattern) − min(pattern))²`

**How the threshold was chosen.** Sort all 24 gate scores (6 methods × 4 clusters) low-to-high. Find the largest gap between consecutive values. Put the threshold at the midpoint of that gap. It's **data-driven, not hand-picked** — the number 22.946 is just where the sorted scores naturally split "methods that land near the pattern" from "methods that don't."

**Why it's universal.** One threshold, applied everywhere. Fixed per-cluster thresholds don't generalize (that was Act 1 of the April 16 notebook — rejected).

### MAE vs MSE (both on raw LT data)
Same pointwise subtraction. MAE uses absolute value, MSE uses the square:
- `MAE = mean(|g_t − p_t|)`
- `MSE = mean((g_t − p_t)²)`

Squaring penalizes big errors disproportionately. For a single gene, MSE ≈ MAE² roughly — rankings are *almost* identical; only genes with one lopsided outlier timepoint flip. Reported best-MAE ≈ 0.0157 and best-MSE ≈ 0.0003 on clusterFour are the same gene on different scales.

---

## Per-cluster story

### clusterOne — spike, moderate magnitude
- 6 normalized methods find spike-ish genes *below* the pattern's level.
- Gated ensemble: 6/6 active (gate didn't silence anyone).
- Raw MAE: genes sit at the right level but are flat — magnitude-matched genes don't actually spike.
- **Verdict**: neither path is clean. No gene has both the right shape and the right level.

### clusterTwo — monotone increase
- Only 0.9% of genes share (+,+,+) — rarest class, shape is actually selective.
- All 6 methods converge on genes that hug both shape and magnitude.
- Gated ensemble: 6/6 active; high consensus (RAB14, SYNGR1, etc.).
- Raw MAE: tight fit, best MAE ≈ 0.03.
- **Verdict**: both paths agree. Easiest cluster.

### clusterThree — dip pattern, tiny magnitude (≈0.04)
- 35.4% of genes share (−,+,−) — the shape-based filter is weakest here.
- 5/6 normalized methods pick genes **10–40× above the pattern** (right shape, wrong level).
- Gated ensemble silences 5, **only sMAPE survives** (1/6 active).
- Raw MAE: genes sit directly on the pattern.
- **Verdict**: magnitude looks closer.

### clusterFour — sharp spike, large magnitude
- 7.4% share (+,−,+).
- 5/6 normalized methods pick **flat genes** — normalization erased the spike.
- Gated ensemble: 6/6 still active. **The gate didn't catch this failure.**
- Raw MAE: captures the spike-shape genes for the first time.
- **Verdict**: magnitude uniquely captures the spike at t=3.

---

## Appendix (raw-space comparison)

Recomputed DTW, Fréchet, MSE, sMAPE on unnormalized LT data to show what the non-scale-invariant methods actually look like without the normalization step.

**Top-20 overlap with raw MAE**:
| Cluster | DTW | Fréchet | MSE | sMAPE |
|---|---|---|---|---|
| clusterOne | 17/20 | 9/20 | 14/20 | 17/20 |
| clusterTwo | 20/20 | 17/20 | 18/20 | 6/20 |
| clusterThree | 20/20 | 13/20 | 16/20 | 17/20 |
| clusterFour | 7/20 | 4/20 | 17/20 | **0/20** |

**Takeaways**:
- **DTW and MSE on raw ≈ raw MAE.** On 4 timepoints, their extra machinery (warping, squared penalty) doesn't distinguish them meaningfully from pointwise MAE.
- **Fréchet is noisier** — max-deviation focus punishes one bad timepoint.
- **sMAPE stands apart.** Percentage-style denominator |g| + |p| breaks down when pattern values approach zero. ClusterFour's near-zero tail (t=6, t=9) makes it pick completely disjoint genes (0/20 overlap).
- **Practical point**: if magnitude-awareness is the goal, raw MAE is the cleanest form of this family.

---

## Anticipated Q&A

**Q: Why 22.946 for the cutoff?**
Not hand-picked. Sort all 24 gate scores, put the threshold at the midpoint of the largest gap — it's where the data naturally splits methods that land near the pattern from methods that don't.

**Q: What's in those 24 scores?**
6 methods × 4 clusters = 24 method-cluster pairs. Each pair gets one relative-centroid-error score.

**Q: Why did clusterFour's gate fail (6/6 active even though methods picked flat genes)?**
Because the relative centroid error is normalized by pattern range squared. ClusterFour's pattern has a large range (tall spike), so `(max − min)²` is big. Even a fairly wrong centroid gets divided down to a moderate error. The gate is fooled by big-range patterns with a few off-genes still within tolerance.

**Q: Is shape similarity or magnitude similarity biologically correct?**
Open question. Shape similarity maps to co-regulation (WGCNA, Eisen tradition) — same upstream signal, magnitudes can differ. Magnitude similarity maps to quantitative/functional similarity — if the constraint represents a signaling threshold or stoichiometric ratio, absolute level matters. Depends on what these constraint patterns actually represent biologically.

**Q: Why only 4 timepoints?**
That's what the experiment produced. With n=4 it's a real limitation — shape classes collapse to 8, and Pearson with n=4 has wide confidence intervals. More timepoints would strengthen shape-based methods.

**Q: Why raw MAE and not raw MSE?**
Near-identical rankings on n=4. MAE is more interpretable (units of expression, not squared). If the group prefers, MSE gives effectively the same picks.

**Q: Why LT (log-transformed) data and not raw counts?**
LT compresses the dynamic range so distance measures aren't dominated by the 2-3 most-expressed genes. It's the convention inherited from Ethan's pipeline. Raw MAE on LT is still "magnitude-aware" — it just compares on a log scale.

**Q: Does sMAPE always win on clusterThree?**
In the gated ensemble, yes — it's the only method that passes the gate. Its top-20 sits on the tiny-magnitude pattern because sMAPE's percentage-style denominator makes it sensitive to relative error, which aligns with the low-magnitude regime. But note sMAPE breaks badly on clusterFour (0/20 overlap with MAE) for the same reason.

**Q: Could we just use raw MAE for everything?**
That's one interpretation of Path B. The catch is clusterOne: no genes at the pattern's magnitude actually follow its spike shape, so raw MAE there picks flat genes at the right level. A hybrid (raw MAE where it works, shape-based for clusterOne) is defensible.

**Q: What's next?**
Depends on the group's answer. If shape → commit to the gated ensemble and maybe refine the gate. If magnitude → raw MAE becomes the method of record, and we'd look at statistical significance and pathway enrichment on the top-N gene lists.

**Q: What did we not do this week?**
- No CSV gene list exports (deliberately scoped out)
- No pathway / GO enrichment (Phase 8)
- No statistical significance testing on n=4
- No clusterFour deep dive (the gate failing on it is a flag for future work)

---

## File map

- `shape_vs_magnitude.ipynb` — main notebook, 6 sections (setup, 4-timepoint problem, Path A, Path B, tradeoffs, literature)
- `raw_space_comparison.ipynb` — recomputes DTW/Fréchet/MSE/sMAPE on unnormalized LT
- `NMF Research Presentation April 23rd MJ.pdf` — 14-slide deck (10 main + 4 appendix)
- `plots/` — all figures for both notebooks
- This file
