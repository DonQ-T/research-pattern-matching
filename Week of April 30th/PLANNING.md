# Week of April 30 — Planning

## Professor feedback (received 2026-04-23)

- **Shape > magnitude is the biological priority.** That resolves the open question from the April 23 presentation.
- **But** raw MAE (Path B) gave the best *shape-wise* visual fits on clusters Two, Three, Four — often better than the shape-based methods.
- **ClusterOne remains the problem child** — no magnitude-matched genes also follow the spike shape.
- **Direction:** try alternative gating metrics (Pearson-based, slope-based, etc.) for the ensemble voting. Goal is **best visual shape fits**, especially on clusterOne.

## Why raw MAE wins shape-wise on 2/3/4

Subtle but important point. The 6 normalized methods end up picking genes with the right sign class but at arbitrary magnitudes — which on low-range clusters (clusterThree) or sharp-spike clusters (clusterFour) means flat genes or wildly scaled ones. Raw MAE picks genes at the right absolute level, and on those clusters the right-level genes *incidentally* also have the right shape.

So: chasing magnitude accidentally chased shape, because the shape-based methods had been neutered by normalization. That's what we want to fix — build a gate that directly optimizes for shape.

## Candidate gating metrics to try

1. **Median Pearson gate** (Apr 16 Gate A, but re-examined now that we know shape is the target)
   - Compute median Pearson correlation between top-20 genes and the pattern.
   - Method passes gate if median > some threshold.
   - Already implemented in April 16 notebook — was rejected for being "too permissive." Worth revisiting: now that we know shape is the goal, "too permissive" might be OK.

2. **Slope-based gate**
   - For each of the 3 transitions (t0→t3, t3→t6, t6→t9), compute the mean slope of the top-20 genes vs. the pattern's slope.
   - Method passes if top-20's mean slopes match the pattern's slopes within tolerance.
   - Captures the *direction and steepness* of change, not just sign.

3. **Sign-class agreement gate**
   - Fraction of top-20 genes that share the constraint pattern's sign class exactly.
   - Simple, interpretable; punishes methods whose picks have wrong-direction transitions.

4. **Hybrid: Pearson + centroid proximity (in normalized space)**
   - Same as Gate C-norm but measure centroid distance *after* normalizing both.
   - Might be redundant with pure Pearson but worth comparing.

5. **Correlation-weighted voting**
   - Instead of silencing methods via a gate, weight each method's votes by how shape-correlated its top-20 are with the pattern.
   - No hard cutoff — smoother.

## Success criterion

Per cluster, which gate + voting scheme produces the top-20 that *visually* best follows the constraint pattern's shape? Especially interested in improvements over raw MAE on **clusterOne**.

Quick comparison table to fill in after experiments:

| Cluster | Raw MAE top-20 | Gate A (Pearson) | Slope gate | Sign-class gate | Correlation-weighted |
|---|---|---|---|---|---|
| clusterOne | ? | ? | ? | ? | ? |
| clusterTwo | ? | ? | ? | ? | ? |
| clusterThree | ? | ? | ? | ? | ? |
| clusterFour | ? | ? | ? | ? | ? |

## Existing assets to reuse

- `../Week of April 9th/metrics_results.pkl` — all 6 methods' scores per cluster
- `../Week of April 23rd/shape_vs_magnitude.ipynb` — Gate C-norm scaffolding, gated ensemble voting function
- `../Week of April 23rd/raw_space_comparison.ipynb` — raw MAE and raw-space metric baselines
- `../Week of April 16th/analysis_cutoffs.ipynb` — original Gate A (median Pearson) code, Gate B (gene spread) for reference

## Out of scope this week

- Pathway / GO enrichment analysis (Phase 8)
- Statistical significance testing on n=4
- Extending to more timepoints (not our data to collect)

## Open questions to revisit

- Should the gate operate on normalized or raw LT data? (Different answers for different candidate gates.)
- Is there a version of the gate that can automatically detect clusterOne's "no good answer" situation and flag it instead of forcing a pick?
