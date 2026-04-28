# Week of April 30: Justification Presentation — Execution Plan

This plan is self-contained. Read this file first, then execute step-by-step. Use TodoWrite to track progress.

## Context

This week's deliverable is to pick the best similarity method **per cluster** (4 clusters total) and export ranked gene lists to two collaborators (Ethan Kyi + Zihan Hei). After running a comparison gallery of 12 PNGs (three y-scale modes × four clusters) in `Week of April 30th/plots/`, Matt reached this visual conclusion:

- **clusterOne**: Pearson and MSE look best (high-magnitude pattern with distinctive up-down-up shape; no genes match magnitude AND shape, so shape is the only recoverable thing per prof's April 23 ruling)
- **clusters Two, Three, Four**: sMAPE and raw-MAE look closest (low-magnitude or zero-touching patterns where rawMAE captures both shape AND magnitude)

Before finalizing the per-cluster picks and shipping the gene lists, Matt wants a presentation + accompanying notes that **convince him both visually and mathematically** that these are the right choices. As a bonus, derive a **dual-criterion gate** that codifies the visual decision process so it's reproducible.

## Project orientation (read these first, in this order)

1. `CLAUDE.md` — project context, methods, data layout, weekly folder convention
2. `C:\Users\mattj\.claude\projects\c--Users-mattj-Desktop-research-pattern-matching\memory\professor_direction.md` — prof's April 23 ruling: shape > magnitude is the priority, but rawMAE happens to give best shape on clusters 2/3/4 because magnitude-matched genes there incidentally have the right shape
3. `Week of April 30th/PLANNING.md` — week-level planning context
4. `Week of April 30th/method_selection.ipynb` — built + executed this week. Read the score-grid outputs in cell `scoregrids` and ensemble outputs in cell `ensembles` for hard numbers.
5. `Week of April 23rd/NOTES.md` — style reference for NOTES.md
6. `Week of April 23rd/shape_vs_magnitude.ipynb` — has helpers worth reusing (sign_pattern, find_gap_threshold, plotting layouts)
7. Past presentations to mimic (Beamer LaTeX style):
   - `Week of April 9th/presentation.tex` + `.pdf`
   - `Week of April 23rd/` — find the .tex / .pdf files
8. The 12 gallery PNGs at `Week of April 30th/plots/gallery_<cluster>_{zoomed,shared,normalized}.png`

## Cluster pattern characteristics (for the math argument)

- **clusterOne**: ~0.85 → 2.55 → 1.85 → 2.05. High magnitude, distinctive up-down-up shape. Sign class shared by ~8.6% of genes.
- **clusterTwo**: ~0 → ~1, monotonic increase. Sign class very rare (~0.9%).
- **clusterThree**: V-shape with t=3 ≈ 0 (a near-zero point that matters). Sign class ~35% of genes.
- **clusterFour**: sharp single spike at t=3, near-flat elsewhere. Sign class ~7.4% of genes.

## Methods (already implemented in `method_selection.ipynb`)

7 base methods (all on log-transformed counts); the first 6 internally min-max normalize before scoring:
1. Pearson correlation
2. DTW
3. Cosine
4. Fréchet
5. MSE (point-wise on normalized)
6. sMAPE
7. raw-MAE (no normalization — operates on raw LT)

Plus 4 gated ensembles (gate metrics: median Pearson, centroid MSE/range², slope match, sign-class agreement).

## Deliverables

- `Week of April 30th/presentation.tex` — Beamer LaTeX matching past weekly style
- `Week of April 30th/presentation.pdf` — compiled, **screenshot-verified for layout** (no overlaps, no overflow)
- `Week of April 30th/NOTES.md` — math + visual rationale + open questions + steelman
- New cell in `method_selection.ipynb` implementing a 5th "dual" gate metric (shape × magnitude)
- New PNGs in `Week of April 30th/plots/`:
  - `dual_scatter_<cluster>.png` (one per cluster — methods plotted in 2D shape×magnitude space)
  - `winner_vs_runnerup_<cluster>.png` (composite showing chosen method's top-20 vs the strongest alternative's top-20, both in zoomed and normalized views)

---

## Step 1 — Style audit of past presentations

Read `Week of April 9th/presentation.tex` and the April 23 deck PDF. Identify:
- Beamer theme + color scheme
- Slide title style + font sizes
- Figure-embedding pattern (filename, caption, sizing)
- How section breaks are introduced
- Typical slide count (~14-16)

Match this style closely. The presentation should look like the next entry in the same series, not a rewrite.

---

## Step 2 — Design and implement the dual-criterion gate

This is the **key technical contribution** of this week's work. The dual gate operationalizes Matt's visual decision process.

### Visual decision process to codify

> "A method is a winner only if it looks good in BOTH the zoomed view (right magnitude) AND the normalized view (right shape). Any method that fails either view falls out of contention."

### Mathematical formulation

For each (method, cluster), compute two axes from the method's top-N=20:

- **Shape axis** — fraction of top-N with `pearsonr(gene_LT, pattern_LT) > 0.9`
  - Range: [0, 1], higher is better
  - Captures BOTH quality (must be high to count) AND consistency (penalizes variance, unlike median Pearson). This is what the eye sees as "clean overlay" vs "fanning" in the normalized view.
- **Magnitude axis** — `mean(|centroid_top_N - pattern|) / pattern.range()` (range-normalized centroid raw-MAE)
  - Range: [0, ∞), lower is better; typically [0, 1] for reasonable picks. Range-normalized so it's comparable across clusters.

### Three combination rules to evaluate (pick the cleanest in the presentation)

1. **Pareto front**: a method is viable iff no other method strictly dominates it on both axes. Among Pareto-optimal, pick the one closest to ideal corner `(shape=1, magnitude=0)` in normalized 2D.
2. **Hard dual gate**: pass iff `shape ≥ S_threshold` AND `magnitude ≤ M_threshold`. Use natural-gap rule (existing notebook helper `find_gap_threshold`) on each axis independently.
3. **Weighted scalar**: `score = 0.7 · shape − 0.3 · normalized_magnitude` (α=0.7 reflects prof's "shape > magnitude" priority).

### Expected behavior (sanity check before writing the math up)

- **clusterOne**: shape methods (Pearson/MSE/Frechet) score high on shape axis, poor on magnitude axis. RawMAE scores well on magnitude axis but mediocre on shape axis (because of fanning visible in normalized view). **No method dominates both** → graceful fallback rule needed: when no method passes both, fall back to highest-shape method (matches prof's priority + Matt's pick).
- **clusters 2/3/4**: rawMAE scores well on BOTH axes (clean dominance). sMAPE may also pass both. → unambiguous dual-gate winner.

If the math gives a different answer than Matt's visual picks, **flag it explicitly** in the notes — either the visual instinct was wrong, or the formula needs adjustment, or there's a defensible alternative worth raising.

### Implementation

Add a new cell to `method_selection.ipynb` after the existing `ensembles` cell. Compute:

- `shape_axis_grid` and `magnitude_axis_grid` — both `pd.DataFrame(index=METHODS_PLUS, columns=cluster_names)`
- `dual_pareto[cluster]` — list of Pareto-optimal methods
- `dual_winner[cluster]` — top pick per cluster under the dual rule

Re-execute the notebook end-to-end to make sure existing outputs still match.

---

## Step 3 — Generate new figures for the presentation

### `dual_scatter_<cluster>.png` (4 files)

For each cluster, a scatter plot in 2D space:
- X-axis: shape score (0-1)
- Y-axis: magnitude score (0-1+, log scale if needed for clusterOne)
- Each method = a labeled point
- Pareto-front line drawn (if using Pareto rule)
- Threshold lines drawn (if using hard gate)
- Winner highlighted with a star marker
- "Ideal corner" at (1, 0) marked with a target symbol

### `winner_vs_runnerup_<cluster>.png` (4 files)

For each cluster, a 2×2 composite:
- Top-left: winner's top-20 in zoomed scale
- Top-right: winner's top-20 in normalized scale
- Bottom-left: runner-up's top-20 in zoomed scale
- Bottom-right: runner-up's top-20 in normalized scale
- Title clearly identifies winner vs runner-up
- Shows visually why winner beats runner-up on both axes (or only on the prof's priority axis for clusterOne)

Choice of runner-up:
- clusterOne: rawMAE (the magnitude-only alternative)
- clusters 2/3/4: Pearson or MSE (the shape-only alternative)

---

## Step 4 — Write the Beamer LaTeX presentation

Target ~14-16 slides. Use the same Beamer theme/color as past weekly decks. Proposed slide structure:

1. **Title** — "Per-Cluster Method Selection: Visual + Mathematical Justification" / Matt Jacob / date / Dr. Guo's group
2. **Recap** — April 23 framing: shape vs magnitude, prof's ruling (shape > magnitude), the paradox (rawMAE wins shape on 2/3/4)
3. **The deliverable** — pick best method per cluster, ship gene CSVs to Ethan + Zihan
4. **Methods compared** — table of 7 base methods + 4 gated ensembles (one row per method, columns: formula, what it minimizes/maximizes, normalization)
5. **Three viewing modes** — explain zoomed / shared / normalized; show the same panel (clusterOne Pearson) in all three modes
6. **Visual scan: clusterOne** — embed `gallery_clusterOne_zoomed.png` and `gallery_clusterOne_normalized.png`. Narrate what jumps out.
7. **Visual scan: clusters 2/3/4** — embed zoomed views for all three clusters. Narrate that rawMAE hugs the line.
8. **The visual decision process** — verbal description of the dual-criterion intuition: "good in zoomed AND good in normalized"
9. **Mathematical formulation** — define shape axis (fraction with Pearson > 0.9) and magnitude axis (centroid MAE / range). Defend each choice in 1-2 sentences.
10. **The dual gate** — formal definition of the chosen combination rule (Pareto / hard gate / weighted) with formula on slide
11. **Per-cluster scatter — clusterOne** — embed `dual_scatter_clusterOne.png`. Show no method passes both gates. Discuss the fallback rule.
12. **Per-cluster scatter — clusters Two/Three/Four** — embed all three. Show rawMAE dominates on each.
13. **Winner vs runner-up — clusterOne** — embed `winner_vs_runnerup_clusterOne.png`. Pearson/MSE shape-perfect; rawMAE chaos in normalized.
14. **Winner vs runner-up — cluster Two/Three/Four** — embed for one representative cluster (probably Three, since it's the V-shape with the zero point that sMAPE handles specially)
15. **Math defense — why these picks are principled** — bullet points connecting each cluster's pattern characteristics to the chosen method's mathematical objective:
    - clusterOne: high-magnitude unrecoverable shape → pure shape methods (Pearson/MSE) are the principled pick because the prof's priority makes them the right fallback when no method passes both gates
    - clusters 2/3/4: rawMAE captures both axes because for these patterns magnitude IS recoverable; sMAPE's symmetric percentage formulation handles near-zero values (clusterThree t=3≈0) gracefully where MSE would be dominated by tiny absolute errors
16. **Steelman** — "what would change my mind":
    - For clusterOne: nothing realistic — already ruled out by prof + visual + math
    - For clusters 2/3/4: if collaborators only care about shape and don't care about magnitude proximity, the "perfect Pearson" picks are equally valid co-regulated genes
    - DTW special note: scored r_med = -0.33 on clusterFour — picked anti-correlated genes. With only 4 timepoints, DTW's warping has too few constraints and can fold the trajectory backwards.
17. **Final picks + deliverable** — table of WINNERS dict per cluster + filename of CSV deliverable
18. *(optional)* **Open questions / next week** — what we'd validate with biology (pathway enrichment), what we'd improve methodologically (more timepoints would help DTW)

Keep slide content terse — Beamer rule of thumb is ≤ 6 bullets per slide, ≤ 10 words per bullet. Move detail to NOTES.md.

---

## Step 5 — Compile + screenshot iteration loop (critical, do not skip)

LaTeX can compile cleanly but produce visually broken output (overlapping text, figures bleeding past margins, slide-overflow). **Verify every page visually before declaring done.**

### Loop

1. Compile: `cd "Week of April 30th" && pdflatex -interaction=nonstopmode presentation.tex` (run twice for cross-references). If `latexmk` is available, prefer `latexmk -pdf presentation.tex`.
2. Render each page to PNG using `pdftoppm`:
   ```
   pdftoppm -png -r 100 presentation.pdf presentation_page
   ```
   This creates `presentation_page-1.png`, `presentation_page-2.png`, etc. at 100 DPI (readable for Read tool).
3. Use the Read tool on each PNG. Specifically check for:
   - **Text overflow** off the right or bottom edge
   - **Figure overflow** past slide boundaries
   - **Overlapping** text and figures
   - **Tiny unreadable** text in tables
   - **Missing figures** (LaTeX placeholder boxes)
   - **Bullet wrap** that breaks into too many lines
4. If any issue: fix in `presentation.tex`, recompile, re-screenshot the affected pages, re-Read.
5. Loop until all pages pass inspection. Typically 2-4 iterations needed.
6. Clean up: delete `presentation_page-*.png` and `*.aux`, `*.log`, `*.nav`, `*.snm`, `*.toc`, `*.out` files at the end. Keep only `.tex` and `.pdf`.

### If `pdftoppm` is unavailable

Fall back to `magick convert -density 100 presentation.pdf presentation_page.png` (ImageMagick). If neither tool is on PATH, ask the user before installing.

### Aspect ratio

Default Beamer is 4:3 (128mm × 96mm). If past presentations used 16:9 (`\documentclass[aspectratio=169]{beamer}`), match that. Check `Week of April 9th/presentation.tex` for the exact preamble.

---

## Step 6 — Write NOTES.md

Style-match `Week of April 23rd/NOTES.md`. Sections:

1. **TL;DR** — 3 bullets summarizing the picks and the dual-gate finding
2. **The visual decision process** — describe what the eye is doing when looking at zoomed vs normalized
3. **The dual-criterion gate** — formal definition + why each axis is principled + how the combination rule was chosen
4. **Why each cluster's pick is right**:
   - clusterOne: walk through the math + visual + prof's priority
   - clusterTwo: dual gate clearly fires for rawMAE
   - clusterThree: special note on sMAPE vs rawMAE for the zero-touching pattern
   - clusterFour: special note on DTW's anti-correlation result + spike rarity
5. **Steelman** — what would change the picks
6. **Open questions for next week** — biology validation, more timepoints, etc.
7. **File pointers** — where the new figures, the updated notebook cell, and the CSVs live

Length target: ~150-300 lines. Should read like a study guide for the deck.

---

## Step 7 — Final WINNERS update + CSV re-export

After the presentation defends the picks, update `WINNERS` in cell `winners` of `method_selection.ipynb`:

```python
WINNERS = {
    'clusterOneLT':   'Pearson',   # or 'MSE' — defend choice in NOTES.md
    'clusterTwoLT':   'rawMAE',
    'clusterThreeLT': 'rawMAE',    # or 'sMAPE' if NOTES.md argues for it
    'clusterFourLT':  'rawMAE',    # or 'sMAPE' if NOTES.md argues for it
}
```

Re-execute cells `winners` through `spotcheck` (sections 6-8) to regenerate CSVs. Verify:
- All 4 sanity checks pass in cell `verify`
- The new CSV top-20 lists match what the presentation argued for
- The combined `final_gene_lists.csv` has the right method names per cluster

Do NOT re-run cells 0-5 — they're slow and the score grids / ensembles / galleries are already correct.

---

## Step 8 — End-state verification

Before declaring done:

1. `presentation.pdf` opens, every page renders cleanly (no overlaps, no overflow) — confirmed via screenshot loop in Step 5
2. `NOTES.md` reads as a coherent study guide
3. `method_selection.ipynb` has the new dual-gate cell and re-runs cleanly end-to-end
4. `deliverables/*.csv` reflect final WINNERS, with the verify cell printing "all checks passed"
5. `plots/` contains the new `dual_scatter_*.png` and `winner_vs_runnerup_*.png` files
6. Git status is sane — no half-written files, no leftover `.aux`/`.log` from LaTeX

Report back to Matt with:
- Path to `presentation.pdf`
- Brief summary of the dual-gate finding (does the math reproduce his visual picks, or did it suggest a tweak?)
- Any open questions raised by the steelman exercise

---

## Constraints

- **Don't modify prior weeks' files** (per CLAUDE.md "preserve Ethan's existing code"). All new work goes in `Week of April 30th/`.
- **Don't create new memories without specific reason** — the existing `professor_direction.md` is sufficient context.
- **Use TodoWrite** to track the 8 steps. Mark each completed as soon as it's done.
- **Don't skip the screenshot loop in Step 5** — visually verifying the PDF is the difference between "compiles" and "presentable."
- **Match past presentation style exactly** — Beamer theme, font sizes, figure sizing, color scheme. This is the 4th deck in a series and should look like the 4th, not a rewrite.
- **If the dual gate's math gives a different answer than Matt's visual picks, raise it** — don't quietly tune the formula to match. Honest disagreement is more useful than fake confirmation.
