# Mathematical Foundations of Similarity Metrics

## Overview

This document explains the mathematical basis for each of the 6 similarity metrics
used in this project, why certain metrics are better suited to gene expression
time-series shape matching, and why some metrics produce overlapping gene rankings.

All metrics operate on **normalized vectors** in $\mathbb{R}^4$ — a pattern vector
$\mathbf{p} = [p_0, p_3, p_6, p_9]$ and a gene expression vector $\mathbf{g} = [g_0, g_3, g_6, g_9]$,
both scaled to $(0, 1]$ via min-max normalization with epsilon offset:

$$x_{\text{norm}} = \frac{x - \min(x) + \varepsilon}{\max(x) - \min(x) + \varepsilon}$$

The $\varepsilon$ (small positive constant) ensures the range is $(0, 1]$ rather than
$[0, 1]$ and prevents division by zero for constant-expression genes.

---

## 1. Pearson Correlation

### Formula

$$r = \frac{\sum_{i}(p_i - \bar{p})(g_i - \bar{g})}{\sqrt{\sum_{i}(p_i - \bar{p})^2 \cdot \sum_{i}(g_i - \bar{g})^2}}$$

where $\bar{p}$ and $\bar{g}$ are the means of $\mathbf{p}$ and $\mathbf{g}$ respectively.

### What It Measures

Pearson correlation is the cosine of the angle between the **mean-centered** versions
of two vectors. By subtracting the mean, it removes any information about the
overall level of expression and isolates the *shape* — the pattern of relative
increases and decreases over time.

The result lies in $[-1, 1]$:
- $r = 1$: identical shape (gene goes up/down at exactly the same timepoints as the pattern)
- $r = 0$: no linear relationship
- $r = -1$: perfectly inverted shape

### Key Properties

- **Scale invariant**: multiplying a gene's expression by any constant does not change r.
  The gene [1, 12, 5, 7] and [100, 1200, 500, 700] yield identical Pearson scores.
- **Shift invariant**: adding a constant to all timepoints does not change r.
  The gene [1, 12, 5, 7] and [101, 112, 105, 107] yield identical scores.
- These invariances mean Pearson captures *shape* regardless of expression magnitude or baseline.

### Why It Suits This Use Case

The biological question is: "does this gene's expression rise and fall at the same
timepoints as the constraint pattern?" That is exactly what mean-centering and
angular comparison measure. A housekeeping gene with constant high expression
[10, 10, 10, 10] has zero variance after mean-centering, yielding r = 0 (or undefined),
correctly indicating no shape match.

### Worked Example

Using the same three genes and pattern from the Cosine section (see Section 2),
all after normalization to $(0, 1]$:

**Pattern:** $\mathbf{p} = [0.0,\; 1.0,\; 0.356,\; 0.479]$ (spike at t=3), mean $\bar{p} = 0.459$

| Gene | Normalized | Mean | Description |
|------|-----------|------|-------------|
| Gene A | $[0.0,\; 1.0,\; 0.371,\; 0.486]$ | 0.464 | Matches pattern shape |
| Gene B | $[1.0,\; 1.0,\; 1.0,\; 1.0]$ | 1.000 | Flat housekeeping gene |
| Gene C | $[0.0,\; 1.0,\; 0.5,\; 0.8]$ | 0.575 | Wobble stretched by normalization |

#### Mean-centering

Subtract each vector's mean to isolate the shape:

- Pattern: $[0.0 - 0.459,\; 1.0 - 0.459,\; 0.356 - 0.459,\; 0.479 - 0.459] = [-0.459,\; 0.541,\; -0.103,\; 0.020]$
- Gene A: $[-0.464,\; 0.536,\; -0.093,\; 0.022]$
- Gene B: $[0,\; 0,\; 0,\; 0]$ (constant gene becomes the zero vector)
- Gene C: $[-0.575,\; 0.425,\; -0.075,\; 0.225]$

#### Computing r for Gene A

$$r = \frac{(-0.459)(-0.464) + (0.541)(0.536) + (-0.103)(-0.093) + (0.020)(0.022)}{\sqrt{0.459^2 + 0.541^2 + 0.103^2 + 0.020^2} \times \sqrt{0.464^2 + 0.536^2 + 0.093^2 + 0.022^2}}$$

$$= \frac{0.213 + 0.290 + 0.010 + 0.0004}{\sqrt{0.506} \times \sqrt{0.510}} = \frac{0.513}{0.711 \times 0.714} = \frac{0.513}{0.508} = \textbf{0.999}$$

Nearly perfect — Gene A's shape almost exactly matches the pattern.

#### Results summary

| Gene | Pearson r | Interpretation |
|------|----------|----------------|
| Gene A | **0.999** | Almost identical shape — best match |
| Gene B | **undefined** | Zero vector after centering, no shape exists |
| Gene C | **0.95** | Relative ups/downs match, but absolute variation was tiny |

**Key takeaway:** Pearson only cares about the *direction* of the centered vector,
not its magnitude. Gene C's original expression barely moved ($[8, 9, 8.5, 8.8]$),
but after normalization and centering, the *pattern* of its tiny wobble matches.
Whether this is biologically meaningful depends on whether a 1-unit fluctuation
in a gene expressed at level 8-9 is real signal or noise.

### Limitations

With $n = 4$ timepoints, the t-test for Pearson significance has $n - 2 = 2$ degrees
of freedom. Under the null hypothesis of no correlation, two random vectors in $\mathbb{R}^4$
have roughly a **20% chance** of $|r| > 0.8$ purely by chance (computed via the
t-distribution: $t = r\sqrt{n-2}/\sqrt{1-r^2}$, $P(|t| > 1.89 \mid df=2) = 0.20$).
This means individual correlation values should not be over-interpreted — the value
lies in the *ranking order* across all ~17,856 genes, not in the absolute magnitude
of $r$ for any single gene.

---

## 2. Cosine Similarity

### Formula

$$\cos(\theta) = \frac{\sum_{i} p_i \cdot g_i}{\sqrt{\sum_{i} p_i^2} \cdot \sqrt{\sum_{i} g_i^2}}$$

### What It Measures

Cosine similarity measures the angle between two vectors in their **original**
(non-centered) form. It is the same dot-product-over-magnitudes calculation as
Pearson, but without first subtracting the mean.

The result lies in $[-1, 1]$ in general, but for non-negative expression data,
it is always in $[0, 1]$.

### Relationship to Pearson

Cosine similarity on raw vectors can be decomposed as:

$$\cos(\mathbf{p}, \mathbf{g}) = \frac{r \cdot \sigma_p \cdot \sigma_g}{\text{rms}(\mathbf{p}) \cdot \text{rms}(\mathbf{g})} + \frac{\bar{p} \cdot \bar{g}}{\text{rms}(\mathbf{p}) \cdot \text{rms}(\mathbf{g})}$$

where $\text{rms}(\mathbf{x}) = \sqrt{\text{mean}(x^2)}$ is the root-mean-square. The first term captures
shape similarity (proportional to Pearson $r$). The second term is a **baseline bias**
that depends on the product of the means.

Pearson correlation equals cosine similarity computed on the mean-centered vectors —
that second term vanishes after centering. When both vectors have large means relative
to their variation (e.g., a high-expression gene with modest fluctuations), the
baseline bias term dominates and cosine returns a high score even if the shapes differ.

### Why It Fails for Shape Matching — Worked Example

Let's compare how Pearson and Cosine evaluate three genes against
clusterOne's constraint pattern. The pattern rises sharply at t=3,
then partially recovers:

**Pattern:** $\mathbf{p} = [1.45,\; 12.0,\; 5.21,\; 6.51]$
(after normalization: $[0.0,\; 1.0,\; 0.356,\; 0.479]$)

| Gene | Expression | Shape | Biology |
|------|-----------|-------|---------|
| Gene A | $[0.5,\; 4.0,\; 1.8,\; 2.2]$ | Spike at t=3, matches pattern | Responsive |
| Gene B | $[10,\; 10,\; 10,\; 10]$ | Flat, constant high expression | Housekeeping |
| Gene C | $[8,\; 9,\; 8.5,\; 8.8]$ | Nearly flat with tiny wobble | Housekeeping |

#### Step 1: Normalization to (0, 1]

- Gene A → $[0.0,\; 1.0,\; 0.371,\; 0.486]$ (looks like the pattern)
- Gene B → $[1.0,\; 1.0,\; 1.0,\; 1.0]$ (all identical — no shape)
- Gene C → $[0.0,\; 1.0,\; 0.5,\; 0.8]$ (looks like it has a spike, but only because min-max stretches a tiny wobble)

#### Step 2: Pearson correlation (mean-centers first)

Pearson subtracts each vector's mean before comparing:

- Pattern centered: $[-0.459,\; 0.541,\; -0.103,\; 0.020]$
- Gene A centered: $[-0.464,\; 0.536,\; -0.093,\; 0.022]$ → **r = 0.999** (nearly identical shape)
- Gene B centered: $[0,\; 0,\; 0,\; 0]$ → **r = undefined** (zero variance, no shape to compare)
- Gene C centered: $[-0.575,\; 0.425,\; -0.075,\; 0.225]$ → **r = 0.95** (the *shape* of the wobble matches, even though in absolute terms it barely moved)

Pearson correctly identifies Gene A as the best match. Gene B is rejected
(no shape). Gene C gets a decent score because its *relative* ups and downs
happen at the same timepoints — whether this is biologically meaningful at
such tiny magnitudes is a separate question.

#### Step 3: Cosine similarity (does NOT mean-center)

Cosine works on the raw normalized vectors:

- Pattern: $[0.0,\; 1.0,\; 0.356,\; 0.479]$, magnitude $= \sqrt{0 + 1 + 0.127 + 0.229} = 1.164$

**Gene A** $[0.0,\; 1.0,\; 0.371,\; 0.486]$:

$$\cos(\theta) = \frac{(0)(0) + (1)(1) + (0.356)(0.371) + (0.479)(0.486)}{1.164 \times 1.175} = \frac{1.365}{1.368} = 0.998$$

**Gene B** $[1.0,\; 1.0,\; 1.0,\; 1.0]$:

$$\cos(\theta) = \frac{(0)(1) + (1)(1) + (0.356)(1) + (0.479)(1)}{1.164 \times 2.0} = \frac{1.835}{2.328} = 0.788$$

**Gene C** $[0.0,\; 1.0,\; 0.5,\; 0.8]$:

$$\cos(\theta) = \frac{(0)(0) + (1)(1) + (0.356)(0.5) + (0.479)(0.8)}{1.164 \times 1.378} = \frac{1.561}{1.604} = 0.973$$

#### The problem

Gene B — the flat housekeeping gene with **zero shape variation** — gets a
cosine similarity of **0.788**. That's high enough to potentially rank it above
thousands of genes that actually have dynamic expression patterns. This happens
because $[1, 1, 1, 1]$ points in a direction that has a positive angle with
almost any vector in the positive quadrant of $\mathbb{R}^4$.

Pearson gives Gene B **undefined/0** — correctly saying "this gene has no shape."
Cosine says **0.788** — incorrectly suggesting it's a decent match.

This is exactly what happens in the real data: Cosine finds housekeeping genes
like EEF1A1 (high, stable expression) in its top-20, while Pearson correctly
excludes them because they have no temporal dynamics matching the pattern.

### When Cosine Is Appropriate

Cosine works well when the baseline level of a signal *is* informative — for
example, in document similarity (TF-IDF vectors), where a zero baseline means
the term is absent. In time-series shape matching, the baseline is not informative,
making cosine a poor choice.

---

## 3. DTW (Dynamic Time Warping)

### Formula

DTW constructs a cost matrix $C$ of size $n \times n$ where:

$$C[i, j] = d(p_i, g_j) + \min\big(C[i{-}1, j],\; C[i, j{-}1],\; C[i{-}1, j{-}1]\big)$$

where $d(p_i, g_j)$ is the pointwise distance (typically squared Euclidean: $(p_i - g_j)^2$).
The DTW distance is $C[n, n]$ — the minimum-cost path through the matrix.

### What It Measures

DTW finds the **optimal temporal alignment** between two time series. Imagine
stretching or compressing the time axis of one series to best align with the other,
subject to the constraint that you can only move forward in time (monotonicity).

A gene that peaks at t = 6 instead of t = 3 would get a high MSE (the peaks
don't line up point-by-point) but potentially a low DTW distance, because DTW
can "warp" t = 6 to align with t = 3.

### Why It Is Largely Redundant at 4 Timepoints

With $n = 4$, the cost matrix is $4 \times 4$. The number of valid warping paths is
extremely limited. While the theoretical maximum shift is 3 positions (mapping
position 1 to position 4), the cost of extreme warps is prohibitive, so practical
warping is constrained to at most 1-2 positions. The warping flexibility that makes
DTW powerful on longer series (15+ points) barely manifests here.

In practice, DTW reduces to a slightly flexible version of pointwise distance,
which is why it shows 16-17/20 overlap with Pearson on clusters One and Three. Both find
genes with the right shape, because there is not enough temporal structure for
DTW's warping to discover genuinely time-shifted genes.

### Worked Example

Same genes and pattern (normalized). We also add a **time-shifted gene** to
show what DTW is designed to catch:

- **Pattern:** $[0.0,\; 1.0,\; 0.356,\; 0.479]$ — spike at t=3
- **Gene A:** $[0.0,\; 1.0,\; 0.371,\; 0.486]$ — spike at t=3 (matches)
- **Gene D:** $[0.0,\; 0.3,\; 1.0,\; 0.4]$ — spike at t=6 (delayed by one timepoint)

Using squared Euclidean distance $d(a,b) = (a - b)^2$:

#### Gene A — aligned match

With Gene A, the optimal path goes straight along the diagonal (no warping needed):

| Step | Pattern point | Gene point | Cost $(p_i - g_i)^2$ |
|------|--------------|------------|----------------------|
| (0,0) | 0.0 | 0.0 | 0.000 |
| (1,1) | 1.0 | 1.0 | 0.000 |
| (2,2) | 0.356 | 0.371 | 0.0002 |
| (3,3) | 0.479 | 0.486 | 0.00005 |

**DTW distance = 0.00025** (essentially zero — perfect alignment)

#### Gene D — delayed spike

Point-by-point (no warping), the errors are large:
- t=3: pattern=1.0 vs gene=0.3 → error = 0.49
- t=6: pattern=0.356 vs gene=1.0 → error = 0.41

But DTW can **warp** pattern's t=3 to align with gene's t=6. The $4 \times 4$ cost
matrix:

|  | g: 0.0 | g: 0.3 | g: 1.0 | g: 0.4 |
|--|--------|--------|--------|--------|
| **p: 0.0** | 0.000 | 0.090 | 1.090 | 1.250 |
| **p: 1.0** | 1.000 | 0.490 | 0.090 | 0.450 |
| **p: 0.356** | 1.127 | 0.493 | 0.505 | 0.092 |
| **p: 0.479** | 1.356 | 0.525 | 0.775 | 0.098 |

Reading the matrix: $C[i,j]$ is the minimum accumulated cost to align
$p_0..p_i$ with $g_0..g_j$. Each cell equals $d(p_i, g_j) + \min(\text{left}, \text{diagonal}, \text{above})$.

**DTW distance = 0.098** — the optimal path warps through (0,0)→(1,1)→(1,2)→(2,3)→(3,3),
effectively aligning pattern's spike (t=3) with gene's spike (t=6).

#### MSE comparison (no warping allowed)

For Gene D point-by-point:

$$\text{MSE} = \frac{(0-0)^2 + (1-0.3)^2 + (0.356-1)^2 + (0.479-0.4)^2}{4} = \frac{0 + 0.49 + 0.41 + 0.006}{4} = 0.227$$

| Gene | DTW | MSE | What happens |
|------|-----|-----|-------------|
| Gene A | 0.00025 | 0.00006 | Both agree: near-perfect match |
| Gene D | **0.098** | **0.227** | DTW finds the shifted spike; MSE penalizes the misalignment |

DTW gives Gene D a score 2.3x better than MSE does (relative to Gene A),
because it can warp the time axis to align the delayed spike. This is DTW's
strength — but with only 4 points, the warping can only shift by 1 position.

### Where DTW Adds Value

On clusterTwo (monotonic increase), DTW diverges from other methods (2-3/20 overlap).
This is the one case where its warping flexibility changes the ranking — it can
identify genes with a slightly delayed increase that point-by-point methods penalize.

With denser time sampling (15+ points), DTW would genuinely capture temporal delays
in gene response (e.g., a gene that activates 2 hours after the stimulus) that
Pearson would miss entirely.

---

## 4. Frechet Distance

### Formula (Discrete)

$$C[i, j] = \max\Big(d(p_i, g_j),\; \min\big(C[i{-}1, j],\; C[i, j{-}1],\; C[i{-}1, j{-}1]\big)\Big)$$

The Fréchet distance is $C[n, n]$.

### What It Measures

The "dog-walking distance": imagine a person walking along curve P and a dog
walking along curve G, connected by a leash. Both can vary their speed but
neither can go backward. The Frechet distance is the length of the shortest
leash that allows both to traverse their respective curves from start to finish.

### Key Difference from DTW

Both DTW and Frechet use dynamic programming with the same recurrence structure.
The critical difference:

- **DTW uses sum** (or equivalently, takes the total accumulated cost)
- **Frechet uses max** (takes the worst single-point deviation)

This means:
- A gene that is slightly off at every timepoint but never terrible:
  **low DTW, moderate Frechet**
- A gene that is perfect at 3 timepoints but way off at 1:
  **low DTW (3 good points dilute the bad one), high Frechet (worst point dominates)**

### Worked Example — Fréchet vs DTW

Two genes illustrate the max-vs-sum difference:

- **Pattern:** $[0.0,\; 1.0,\; 0.356,\; 0.479]$
- **Gene E:** $[0.05,\; 0.90,\; 0.40,\; 0.55]$ — slightly off at every timepoint, but never far
- **Gene F:** $[0.0,\; 1.0,\; 0.356,\; 0.9]$ — perfect at 3 timepoints, but way off at t=9

#### Pointwise distances $(p_i - g_i)^2$

| Timepoint | Gene E error | Gene F error |
|-----------|-------------|-------------|
| t=0 | $(0.0 - 0.05)^2 = 0.0025$ | $(0.0 - 0.0)^2 = 0.0$ |
| t=3 | $(1.0 - 0.90)^2 = 0.01$ | $(1.0 - 1.0)^2 = 0.0$ |
| t=6 | $(0.356 - 0.40)^2 = 0.002$ | $(0.356 - 0.356)^2 = 0.0$ |
| t=9 | $(0.479 - 0.55)^2 = 0.005$ | $(0.479 - 0.9)^2 = 0.177$ |

#### DTW (sum of costs along diagonal path)

Since both genes align best without warping, DTW takes the diagonal path:

- **Gene E DTW:** $0.0025 + 0.01 + 0.002 + 0.005 = 0.020$
- **Gene F DTW:** $0.0 + 0.0 + 0.0 + 0.177 = 0.177$

DTW says Gene E is **8.9x closer** — the three perfect points in Gene F don't
compensate for the one bad point, because the errors still accumulate.

#### Fréchet (max cost along the path)

Fréchet takes only the worst single deviation:

- **Gene E Fréchet:** $\max(0.0025, 0.01, 0.002, 0.005) = 0.01$
- **Gene F Fréchet:** $\max(0.0, 0.0, 0.0, 0.177) = 0.177$

Fréchet says Gene E is **17.7x closer** — it punishes Gene F's single bad
timepoint even more harshly than DTW does, because it *only* looks at the
worst point.

#### Summary

| Metric | Gene E | Gene F | Which is harsher on Gene F? |
|--------|--------|--------|----------------------------|
| DTW | 0.020 | 0.177 | 8.9x worse |
| Fréchet | 0.01 | 0.177 | **17.7x worse** |
| MSE | 0.005 | 0.044 | 8.9x worse |

Fréchet is the most conservative metric — a gene must be close at **every**
timepoint to score well. DTW and MSE are more forgiving because they average
out the errors.

### Behavior on This Data

With only 4 timepoints, the max-vs-sum distinction creates minimal divergence.
The largest error term tends to dominate both metrics. This is reflected in the
data: Fréchet and MSE show 19-20/20 overlap on clusters Three and Four, and
15-20/20 overlap across other clusters.

---

## 5. MSE (Mean Squared Error)

### Formula

$$\text{MSE} = \frac{1}{n} \sum_{i}(p_i - g_i)^2$$

### What It Measures

The simplest possible distance metric: the average squared difference between
corresponding timepoints. No alignment, no normalization tricks, no angular
measurement. Just "how far apart are these two vectors, point by point?"

MSE is proportional to the square of the Euclidean distance:

$$\text{MSE} = \frac{1}{n} \|\mathbf{p} - \mathbf{g}\|^2$$

### Worked Example

Same genes, same pattern:

- **Pattern:** $[0.0,\; 1.0,\; 0.356,\; 0.479]$
- **Gene A:** $[0.0,\; 1.0,\; 0.371,\; 0.486]$ (matches shape)
- **Gene B:** $[1.0,\; 1.0,\; 1.0,\; 1.0]$ (flat housekeeping)
- **Gene C:** $[0.0,\; 1.0,\; 0.5,\; 0.8]$ (wobble)

#### Computing MSE for each gene

**Gene A:**

$$\text{MSE} = \frac{(0.0-0.0)^2 + (1.0-1.0)^2 + (0.356-0.371)^2 + (0.479-0.486)^2}{4} = \frac{0 + 0 + 0.000225 + 0.000049}{4} = 0.00007$$

**Gene B:**

$$\text{MSE} = \frac{(0.0-1.0)^2 + (1.0-1.0)^2 + (0.356-1.0)^2 + (0.479-1.0)^2}{4} = \frac{1.0 + 0 + 0.415 + 0.271}{4} = 0.422$$

**Gene C:**

$$\text{MSE} = \frac{(0.0-0.0)^2 + (1.0-1.0)^2 + (0.356-0.5)^2 + (0.479-0.8)^2}{4} = \frac{0 + 0 + 0.0208 + 0.103}{4} = 0.031$$

#### Results compared to Pearson

| Gene | MSE | Pearson r | MSE rank | Pearson rank |
|------|-----|----------|----------|-------------|
| Gene A | **0.00007** | **0.999** | 1st (best) | 1st (best) |
| Gene C | 0.031 | 0.95 | 2nd | 2nd |
| Gene B | 0.422 | undefined | 3rd (worst) | 3rd (worst) |

MSE and Pearson agree on the ranking here. But notice what drives the scores:

- **Pearson** gives Gene C a high score (0.95) because the *relative shape* matches.
  It doesn't care that Gene C's actual values $[0.0, 1.0, 0.5, 0.8]$ deviate from
  the pattern $[0.0, 1.0, 0.356, 0.479]$ by a fair amount at t=6 and t=9.
- **MSE** penalizes Gene C for those absolute deviations (0.031 vs 0.00007 for Gene A —
  a 440x difference), even though the shape is similar.

This is where the 2-4/20 disagreements between Pearson and MSE come from: genes
with the right shape but slightly different proportions get high Pearson but
mediocre MSE.

### Sensitivity to Scale

Unlike Pearson, MSE is sensitive to both shape *and* scale. Without normalization,
MSE would rank genes primarily by expression level — high-expression genes
would always have large squared differences. This is why the $(0, 1]$ normalization
is essential: it puts all genes on a comparable scale, allowing MSE to function
as a shape metric.

### Relationship to Pearson

The full MSE-Pearson decomposition is:

$$\text{MSE} = \text{Var}(\mathbf{p}) + \text{Var}(\mathbf{g}) - 2r\,\sigma_p\sigma_g + (\bar{p} - \bar{g})^2$$

When the two vectors have **equal variance** AND **equal mean**, this simplifies to:

$$\text{MSE} = 2\sigma^2(1 - r)$$

In practice, min-max normalized vectors do not have identical means or variances,
so the full decomposition applies. The $(\bar{p} - \bar{g})^2$ term means two vectors
can have high Pearson $r$ but still differ in MSE if their means differ — this
explains some of the 2-4/20 disagreements between the methods. Nevertheless,
the $r$ term dominates for most genes, which is why Pearson and MSE show 16-18/20
overlap in the data. They are measuring nearly the same thing from complementary
perspectives — Pearson via angular similarity, MSE via Euclidean distance.

### Value as Cross-Validation

Despite the overlap, MSE provides independent confirmation from a different
mathematical framework. If Pearson and MSE agree on the top genes, confidence
is higher than from either method alone. The few genes where they disagree
(typically 2-4 out of 20) can reveal cases where the relationship is more
nuanced than simple linear correlation.

---

## 6. sMAPE (Symmetric Mean Absolute Percentage Error)

### Formula

$$\text{sMAPE} = \frac{1}{n} \sum_{i} \frac{2\,|p_i - g_i|}{|p_i| + |g_i|}$$

Each term lies in $[0, 2]$. The "symmetric" form divides by the average of both
values rather than just the reference, avoiding the asymmetry of standard MAPE.

### What It Measures

sMAPE measures **relative** (percentage) error at each timepoint. A difference
of 0.01 when both values are around 0.01 counts the same as a difference of 1.0
when both values are around 1.0. This makes it scale-sensitive at each individual
timepoint, weighting errors proportionally to the local magnitude.

### Worked Example

Using the same genes against the clusterOne pattern, plus clusterThree to show
the zero-value failure:

**clusterOne pattern (normalized):** $\mathbf{p} = [0.0,\; 1.0,\; 0.356,\; 0.479]$

#### sMAPE for each gene against clusterOne

**Gene A** $[0.0,\; 1.0,\; 0.371,\; 0.486]$:

| Timepoint | $\|p_i - g_i\|$ | $\|p_i\| + \|g_i\|$ | sMAPE term |
|-----------|----------|-------------|-----------|
| t=0 | 0.0 | 0.0 + 0.0 = 0.0 | 0/0 → 0 (both zero) |
| t=3 | 0.0 | 1.0 + 1.0 = 2.0 | 0.0 |
| t=6 | 0.015 | 0.356 + 0.371 = 0.727 | $2(0.015)/0.727 = 0.041$ |
| t=9 | 0.007 | 0.479 + 0.486 = 0.965 | $2(0.007)/0.965 = 0.015$ |

$$\text{sMAPE} = \frac{0 + 0 + 0.041 + 0.015}{4} = \textbf{0.014}$$

**Gene B** $[1.0,\; 1.0,\; 1.0,\; 1.0]$:

| Timepoint | $\|p_i - g_i\|$ | $\|p_i\| + \|g_i\|$ | sMAPE term |
|-----------|----------|-------------|-----------|
| t=0 | 1.0 | 0.0 + 1.0 = 1.0 | $2(1.0)/1.0 = 2.0$ |
| t=3 | 0.0 | 1.0 + 1.0 = 2.0 | 0.0 |
| t=6 | 0.644 | 0.356 + 1.0 = 1.356 | $2(0.644)/1.356 = 0.950$ |
| t=9 | 0.521 | 0.479 + 1.0 = 1.479 | $2(0.521)/1.479 = 0.704$ |

$$\text{sMAPE} = \frac{2.0 + 0 + 0.950 + 0.704}{4} = \textbf{0.914}$$

**Gene C** $[0.0,\; 1.0,\; 0.5,\; 0.8]$:

| Timepoint | $\|p_i - g_i\|$ | $\|p_i\| + \|g_i\|$ | sMAPE term |
|-----------|----------|-------------|-----------|
| t=0 | 0.0 | 0.0 | 0/0 → 0 |
| t=3 | 0.0 | 2.0 | 0.0 |
| t=6 | 0.144 | 0.856 | $2(0.144)/0.856 = 0.336$ |
| t=9 | 0.321 | 1.279 | $2(0.321)/1.279 = 0.502$ |

$$\text{sMAPE} = \frac{0 + 0 + 0.336 + 0.502}{4} = \textbf{0.210}$$

#### Notice the problem at t=0

When $p_i = 0$ and $g_i = 0$ (Gene A and Gene C at t=0), we get $0/0$. When
$p_i = 0$ but $g_i = 1.0$ (Gene B at t=0), the sMAPE term maxes out at 2.0 —
the largest possible value for a single timepoint.

In the clusterOne example above, this doesn't cause a huge problem because the
other timepoints also produce sizable sMAPE terms. But on **clusterThree** — where
the pattern is $[\sim1, \sim0, \sim0.91, \sim0.87]$ — the zero at t=3 contributes a
fixed 2.0 for *every gene*, while the other three near-constant timepoints produce
only tiny sMAPE differences between genes. The ranking ends up decided by noise
at the non-informative timepoints, because the one biologically important timepoint
(the dip) has been reduced to a constant. That's the real failure case, detailed below.

#### All methods compared

| Gene | Pearson | MSE | DTW | Cosine | sMAPE | Description |
|------|---------|-----|-----|--------|-------|-------------|
| Gene A | 0.999 | 0.00007 | 0.00025 | 0.998 | 0.014 | Shape match ✓ |
| Gene B | undef | 0.422 | 1.686 | 0.788 | 0.914 | Flat ✗ |
| Gene C | 0.95 | 0.031 | 0.124 | 0.973 | 0.210 | Wobble ≈ |

All six methods agree: Gene A is the best match, Gene B is the worst. But they
disagree on *how bad* Gene B is and *how good* Gene C is — and those disagreements,
scaled across 17,856 genes, are what produce the different top-20 lists.

### Why It Fails on clusterThree

clusterThree's pattern has a zero at t = 3 (after normalization: [~1, ~0, ~0.91, ~0.87]).

When $p_i = 0$ (or very near zero):

$$\text{sMAPE}_i = \frac{2\,|0 - g_i|}{|0| + |g_i|} = \frac{2\,|g_i|}{|g_i|} = 2.0$$

This equals exactly 2.0 regardless of what $g_i$ is (as long as $g_i \neq 0$). When both
$p_i = 0$ and $g_i = 0$, the formula becomes $0/0$ (undefined); implementations typically
return 0 for this case. Either way, the most informative timepoint — the dip that
defines clusterThree's shape — contributes a constant (or near-constant) to every
gene's score and provides zero discriminating power.

sMAPE then ranks genes based solely on percentage error at the other 3 near-constant
timepoints, completely missing the biologically important feature. This is why sMAPE
shows **0 overlap** with every other method on clusterThree.

### When sMAPE Is Appropriate

sMAPE is designed for forecasting evaluation, where all values are well away from
zero and proportional accuracy matters more than absolute accuracy. It is not
appropriate for time-series shape matching when the reference series passes
through or near zero.

---

## Why Methods Cluster Together

The overlap matrices from the analysis reveal a clear grouping structure:

### Group 1: Shape/Distance Metrics (Pearson, MSE, Frechet)

After normalization, all three measure how close a gene's trajectory is to the
pattern's trajectory. Pearson does it via angular comparison (after mean-centering),
MSE via squared Euclidean distance, and Frechet via worst-case deviation. Because
the approximate relationship $\text{MSE} \approx 2\sigma^2(1-r)$ holds for normalized data, and
Fréchet is dominated by the same large-error terms as MSE at $n = 4$, these three
consistently show 15-20/20 overlap.

### Group 2: Alignment Metric (DTW)

DTW overlaps with Group 1 when temporal warping is irrelevant (clusterThree:
17-19/20 overlap). It diverges when warping matters (clusterTwo: 2-3/20 overlap).
It measures something genuinely different (time-flexible distance), but with 4
timepoints that difference rarely manifests.

### Group 3: Context-Dependent (Cosine)

Cosine's behavior varies significantly across clusters because its magnitude bias
interacts differently with each pattern's geometry:
- **clusterThree** (near-flat pattern): Cosine-Pearson overlap = **20/20** (perfect).
  When the pattern has low variance, mean-centering has less effect, so Cosine and
  Pearson converge.
- **clusterOne** (spike pattern): Cosine-Pearson overlap = **2/20**. The spike creates
  high variance, so mean-centering matters and Cosine finds different (less relevant) genes.
- **clusterFour** (extreme spike): Cosine clusters with Frechet/MSE (19-20/20 overlap)
  rather than with Pearson (0/20). See "The clusterFour Problem" below.

This inconsistency makes Cosine unreliable as a primary metric — its rankings depend
on data geometry in ways that don't track biological relevance.

### Group 4: Percentage-Based (sMAPE)

Completely isolated when zeros are present (clusterThree: 0 overlap with everything).
Moderate overlap in other clusters (clusterOne: 8-9/20), but measuring something
sufficiently different to not cluster with any other group.

---

## The clusterFour Problem

clusterFour's pattern is [0.009, 0.638, 0.018, 0.024] — a 70x spike at t = 3.
After normalization, this is effectively [0, 1, 0, 0]: a near-delta function.

Matching a delta function requires exact agreement at one timepoint (t = 3), while
the other 3 timepoints are all near zero and contribute almost no discriminating
information. Each method projects "closeness to a spike" differently:

- **Pearson**: any gene with a relative peak at t = 3, regardless of magnitude (scale invariant)
- **MSE/Frechet**: dominated by the squared/max error at the t = 3 value
- **Cosine**: dominated by the t = 3 component of the vector
- **DTW**: can warp neighboring timepoints toward t = 3, changing which genes win
- **sMAPE**: near-zero values at t = 0, 6, 9 cause percentage calculations to blow up

The overlap matrix reveals a surprising structure: while Pearson, DTW, and sMAPE
are completely isolated from each other (0/20 overlap), **Cosine, Frechet, and MSE
form a tight cluster** (19-20/20 overlap). This happens because all three are
dominated by the magnitude of each gene's value at t = 3: Cosine by the dot product
component, Frechet by the max deviation at the spike, and MSE by the squared error
at the spike. Since the spike is so extreme, the ranking effectively reduces to
"which genes have the largest expression at t = 3" — and all three metrics agree
on that ordering.

Pearson diverges because it is scale-invariant: it finds genes with a *relative*
peak at t = 3 regardless of absolute magnitude, producing a different ranking.
DTW diverges because its warping flexibility changes which neighbors get aligned
to the spike.

**Resolution**: this cluster needs denser time sampling (e.g., t = 0, 1, 2, 3, 4, 5, 6, 9)
to provide enough structure around the spike for all methods to converge. The
Cosine/Frechet/MSE agreement here is based on a single-timepoint ranking, not
genuine shape matching.

---

## Recommendations for This Use Case

1. **Primary metric: Pearson correlation** — most interpretable, directly captures
   shape similarity, and its invariance properties are exactly what biological
   pattern matching requires.

2. **Secondary validation: MSE or Frechet** — provides independent distance-based
   confirmation. When Pearson and MSE agree, confidence is high.

3. **Skip for this data: Cosine** (magnitude bias makes it find housekeeping genes
   instead of shape-matching genes) and **sMAPE** (percentage-based formulation
   breaks down at zero or near-zero values).

4. **DTW**: marginal benefit at 4 timepoints. Worth including in the ensemble for
   completeness, and would add substantial value with denser time series.

5. **Ensemble strategy**: Frechet and MSE are near-redundant (19-20/20 overlap on
   most clusters), so a 4-method ensemble (Pearson + DTW + Frechet + MSE) with
   a >= 3/4 vote threshold effectively requires agreement between Pearson (or DTW)
   and the Frechet/MSE pair. A 3-method ensemble (Pearson + DTW + MSE) with >= 2/3
   agreement may be more honest about the effective degrees of freedom. Either way,
   genes that appear in multiple independent methods are the highest-confidence candidates.

6. **clusterFour**: needs either more timepoints or a specialized peak-detection
   approach rather than whole-series similarity matching.
