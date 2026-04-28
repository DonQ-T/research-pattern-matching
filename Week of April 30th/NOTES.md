# Week of April 30 — Study Notes

Study guide for the rawMAE-shape-test framework and the per-cluster picks.

---

## TL;DR

1. **Picks**: clusterOne → **Pearson**; clusterTwo, Three, Four → **rawMAE**.
2. **The rule**: one question per cluster — *Is the magnitude recoverable?* — answered by one number, **rawMAE's top-20 shape score** (fraction with `Pearson(gene, pattern) > 0.9`).
   - High (≥0.9) → use rawMAE; it wins shape AND magnitude.
   - Low (<0.9) → fall back to **Pearson**, the principled scale-invariant shape metric.
3. **The numbers**: `0.35, 1.00, 1.00, 1.00`. clusterOne is the lone "no" — magnitude unrecoverable, Pearson takes over.

---

## The framework

### Why this rule

April 23 settled the priority: shape > magnitude. The April 23 paradox was that rawMAE — a *magnitude* metric — was outperforming the explicit shape methods on three of four clusters. Why? Because the shape methods ran on normalized data, which erased magnitude. They'd find shape-correct genes at the wrong amplitude. rawMAE found amplitude-correct genes that *incidentally* had the right shape too, so it won both axes by accident.

The rawMAE-shape-test directly checks whether the accident is happening on a given cluster:
- If rawMAE's top-20 are also shape-correct (Pearson > 0.9 vs pattern), the accident is real → use rawMAE.
- If rawMAE's top-20 are NOT shape-correct, the accident isn't there → fall back to a method whose math actually targets shape (Pearson).

### Why Pearson as the fallback

When we fall back, we need the metric whose mathematical objective is shape. That's Pearson:

```
r(g, p) = cov(g, p) / (σ_g · σ_p)
        = Σ (g_t − ḡ)(p_t − p̄) / sqrt(Σ(g_t − ḡ)² · Σ(p_t − p̄)²)
```

Pearson explicitly mean-centers and divides by std. Its score IS scale invariance. The other "shape" methods are misnamed:

| Method | What it really optimizes | "Shape" only because... |
|---|---|---|
| MSE | Point-wise squared distance | data was normalized first |
| Fréchet | Max curve deviation | data was normalized first |
| sMAPE | Symmetric percentage error | data was normalized first |
| DTW | Time-warped distance | data was normalized first; n=4 lets it fold trajectories |
| Cosine | Vector angle from origin | similar to Pearson but does NOT mean-center |

If we strip the upstream normalization, all of those become magnitude metrics. Pearson stays Pearson — its math doesn't depend on a preprocessing step.

### Why this story is simpler than the dual gate

Last draft used a 2D Pareto rule on (shape, magnitude) axes. It was technically correct but hid the simple truth: **rawMAE works when its picks are shape-correct, and Pearson is the principled fallback when they aren't.** One question, one number, four answers. The dual-gate cell stays in the notebook as supporting analysis (it gives the same answer up to a clusterOne disagreement that we now interpret correctly), but the deck and shipping rule lean on the simpler test.

---

## Per-cluster picks

### clusterOne — pattern: 0.85 → 2.55 → 1.85 → 2.05  (high magnitude, up-down-up)

**Pick: Pearson.** rawMAE shape = **0.35** — only 7 of rawMAE's top-20 are shape-correct. The "free shape" trick fails here because no genes match both magnitude AND the up-down-up shape.

**Why Pearson over the other shape-on-normalized methods:**

1. **Mathematical principle.** Pearson is the only metric whose math explicitly defines shape correlation. The others are distance/angle metrics that became shape-like only because of upstream normalization. Without that step they're magnitude metrics.
2. **Empirical consensus.** clusterOne top-20 overlap with Pearson:
   - Fréchet: 18/20  •  MSE: 17/20  •  DTW: 13/20  •  sMAPE: 4/20  •  Cosine: 0/20
   Pearson, Fréchet, MSE form a tight consensus (different math, same gene picks: NACA, AC008875.3, MDFIC, RTRAF, FAM89B). Pearson is the *principled* member of that consensus.
3. **Cosine and sMAPE pick fundamentally different genes.** They're shape-correct in their own way, but their math doesn't directly target shape. No reason to prefer them.
4. **Gated ensembles aggregate these methods.** They reproduce the Pearson/Fréchet/MSE consensus; they don't improve on it.

**The mean-shift visualization:** translate each gene additively so its mean equals the pattern's mean. This removes the magnitude excuse and forces a pure-shape judgment. Pearson's top-20 mean-shifted trace the up-down-up shape beautifully. rawMAE's top-20 mean-shifted are still chaos — wrong shape isn't fixable by translation. See `clusterOne_meanshift_compare.png`.

### clusterTwo — pattern: 0.37 → 0.80 → 1.22 → 1.78  (steady increase)

**Pick: rawMAE.** rawMAE shape = **1.00** — every one of rawMAE's top-20 has Pearson > 0.9 vs the pattern. Magnitude AND shape both nailed.

Pearson alternative on this cluster has the same shape but spreads ~5× further from the pattern in absolute LT distance. rawMAE picks tighter, magnitude-aligned co-regulated genes.

### clusterThree — pattern: 0.041 → 0 → 0.037 → 0.035  (V-shape, t=3 ≈ 0)

**Pick: rawMAE.** rawMAE shape = **1.00** — same story. The V's tiny magnitude (range ≈ 0.04) means normalized methods' shape-correct picks live at LT 0.5–1.0 (5× the pattern's range). rawMAE finds genes literally on the V.

### clusterFour — pattern: 0.009 → 0.638 → 0.018 → 0.024  (sharp spike at t=3)

**Pick: rawMAE.** rawMAE shape = **1.00**.

**DTW failure note**: DTW scored r_med = −0.33 here — it picked anti-correlated genes. With n=4 timepoints, DTW's warping has too few constraints and can fold the trajectory backwards. **Don't use DTW on n=4 spike patterns.**

---

## Steelman — what would change the picks

| Cluster | Current pick | What would flip it |
|---|---|---|
| clusterOne | Pearson | Fréchet or MSE: same 17–18/20 picks, would barely change the deliverable. Use Pearson for the principled-by-math reason. |
| clusterOne | Pearson | Cosine or sMAPE: different gene picks; would need a biological reason to prefer them. None offered today. |
| clusters 2–4 | rawMAE | Pearson if collaborators only care about shape (not magnitude). Different biological question, different answer. |

---

## Open questions for next week

- **Pathway / GO enrichment** on the four shipped top-20 lists — does each cluster's gene set hit a coherent biological process? Missing biological validation.
- **DTW on more timepoints**: the n=4 anti-correlation failure on clusterFour is a warping over-fit. Worth retrying once denser temporal sampling exists.
- **Generalize the rawMAE-shape-test** to a second NMF run with different constraint patterns. If it picks the right method automatically, the rule is shippable as a default selector beyond TCR.

---

## File map

- `presentation.tex` / `.pdf` — 23-slide Beamer deck
- `method_selection.ipynb` — section 4.5 has the dual-gate analysis (supporting); the `winners` cell prints the rawMAE-shape-test result
- **New plots**:
  - `clusterOne_meanshift_compare.png` — Pearson vs rawMAE in zoomed and mean-shifted views (the "translate up" visual)
  - `clusterOne_shape_methods_meanshift.png` — 6 shape methods' top-20 mean-shifted, with overlap counts
  - `clusterOne_pearson_vs_each_meanshift.png` — paired Pearson vs each alternative
  - `winner_vs_alt_<cluster>.png` (×4) — winner vs runner-up, zoomed and normalized, aligned to the new picks
- `plots/gallery_<cluster>_{zoomed,shared,normalized}.png` — 12 gallery files
- `deliverables/<cluster>_ranking.csv` — full 17,856-gene ranking under each cluster's winner
- `deliverables/final_gene_lists.csv` — combined; 71,424 rows
- `PLANNING.md` — week-level planning context
- `PRESENTATION_PLAN.md` — original execution plan (kept for history; superseded by the simpler framework above)
