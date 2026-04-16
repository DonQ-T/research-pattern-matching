"""Per-gene breakdown figures: show the per-timepoint contributions
that aggregate to each final metric score.

Layout per metric: 3 columns (one per gene) x 2 rows
  row 0 — pattern + gene curves overlay
  row 1 — bar chart of per-timepoint contribution (the inputs)
  bottom text — the aggregation step + final score
"""
import numpy as np
import matplotlib.pyplot as plt

x = [0, 3, 6, 9]

# ----------------- Standard test vectors (normalized) -----------------
pattern = np.array([0.0, 1.0, 0.356, 0.479])

gene_A = np.array([0.0, 1.0, 0.371, 0.486])  # shape match
gene_B = np.array([1.0, 1.0, 1.0, 1.0])      # housekeeping
gene_C = np.array([0.0, 1.0, 0.5, 0.8])      # wobble
gene_D = np.array([0.0, 0.3, 1.0, 0.4])      # delayed spike
gene_E = np.array([0.05, 0.90, 0.40, 0.55])  # small errors everywhere
gene_F = np.array([0.0, 1.0, 0.356, 0.9])    # 3 perfect, 1 disaster


# ----------------- Reusable plotters -----------------
def plot_curves(ax, p, g, gene_label, title):
    ax.plot(x, p, 'o-', color='crimson', lw=2.2, label='Pattern', ms=7)
    ax.plot(x, g, 's-', color='steelblue', lw=2.2, label=gene_label, ms=7)
    ax.set_xticks(x)
    ax.set_xlabel('Timepoint (h)')
    ax.set_ylim(-0.1, 1.2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title(title, fontsize=10, fontweight='bold')


def plot_contrib_bars(ax, positions, contribs, xticklabels, formula_text,
                      ylabel, summary_text, max_idx=None):
    pos_color = 'steelblue'
    neg_color = 'lightcoral'
    colors = [pos_color if c >= 0 else neg_color for c in contribs]
    if max_idx is not None:
        colors[max_idx] = 'red'
    ax.bar(positions, contribs, color=colors, edgecolor='black', linewidth=0.7)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xticks(positions)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('Timepoint (h)')
    ax.set_ylabel(ylabel)
    ax.set_title(formula_text, fontsize=9)
    ax.grid(alpha=0.3, axis='y')
    for xi, c in zip(positions, contribs):
        ax.text(xi, c, f'{c:.4f}', ha='center',
                va='bottom' if c >= 0 else 'top', fontsize=8, fontweight='bold')
    ax.text(0.5, -0.40, summary_text, transform=ax.transAxes,
            fontsize=9.5, ha='center', va='top', color='darkred',
            fontweight='bold',
            bbox=dict(facecolor='lightyellow', edgecolor='gray', boxstyle='round,pad=0.4'))


def finalize(fig, suptitle, save):
    fig.suptitle(suptitle, fontsize=12.5, fontweight='bold')
    plt.tight_layout()
    plt.subplots_adjust(top=0.92, hspace=0.55, bottom=0.16)
    plt.savefig(f'plots/{save}.png', dpi=130, bbox_inches='tight')
    plt.close()


# =====================================================================
# 1. PEARSON  -- per-timepoint cross-product of mean-centered values
# =====================================================================
def pearson_breakdown(p, g):
    p_c = p - p.mean()
    g_c = g - g.mean()
    cross = p_c * g_c
    num = float(cross.sum())
    den = float(np.sqrt((p_c ** 2).sum()) * np.sqrt((g_c ** 2).sum()))
    r = num / den if den > 0 else float('nan')
    return cross, num, den, r


genes_pearson = [('Gene A (match)', gene_A),
                 ('Gene C (wobble)', gene_C),
                 ('Gene B (flat)', gene_B)]
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
for i, (name, g) in enumerate(genes_pearson):
    cross, num, den, r = pearson_breakdown(pattern, g)
    plot_curves(axes[0, i], pattern, g, name, name)
    formula = r'contribution at $t_i$:  $(p_i - \bar p)(g_i - \bar g)$'
    if np.isnan(r):
        summary = (f'numerator (sum) = {num:.4f}\n'
                   f'denominator = {den:.4f}\n'
                   f'r = 0 / 0  =  undefined')
    else:
        summary = (f'numerator (sum) = {num:.4f}\n'
                   f'denominator (mag $\\times$ mag) = {den:.4f}\n'
                   f'r = {num:.4f} / {den:.4f}  =  {r:.4f}')
    plot_contrib_bars(axes[1, i], x, cross, [str(t) for t in x],
                      formula, 'cross-product', summary)
finalize(fig,
         'Pearson breakdown — sum the cross-products of mean-centered values, then normalize',
         'method_pearson_breakdown')


# =====================================================================
# 2. COSINE  -- per-timepoint product, divided by vector magnitudes
# =====================================================================
def cosine_breakdown(p, g):
    prods = p * g
    num = float(prods.sum())
    den = float(np.linalg.norm(p) * np.linalg.norm(g))
    cos = num / den if den > 0 else float('nan')
    return prods, num, den, cos


genes_cosine = [('Gene A (match)', gene_A),
                ('Gene C (wobble)', gene_C),
                ('Gene B (flat)', gene_B)]
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
for i, (name, g) in enumerate(genes_cosine):
    prods, num, den, cos = cosine_breakdown(pattern, g)
    plot_curves(axes[0, i], pattern, g, name, name)
    formula = r'contribution at $t_i$:  $p_i \cdot g_i$'
    summary = (f'numerator (sum of products) = {num:.4f}\n'
               f'denominator $\\|p\\| \\cdot \\|g\\|$ = {den:.4f}\n'
               f'cos = {num:.4f} / {den:.4f}  =  {cos:.4f}')
    plot_contrib_bars(axes[1, i], x, prods, [str(t) for t in x],
                      formula, 'pointwise product', summary)
finalize(fig,
         'Cosine breakdown — sum the raw pointwise products, divide by the magnitudes',
         'method_cosine_breakdown')


# =====================================================================
# 5. MSE  -- per-timepoint squared error
# =====================================================================
def mse_breakdown(p, g):
    sq = (p - g) ** 2
    s = float(sq.sum())
    return sq, s, s / len(sq)


genes_mse = [('Gene A (match)', gene_A),
             ('Gene C (wobble)', gene_C),
             ('Gene B (flat)', gene_B)]
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
for i, (name, g) in enumerate(genes_mse):
    sq, s, mse = mse_breakdown(pattern, g)
    plot_curves(axes[0, i], pattern, g, name, name)
    formula = r'contribution at $t_i$:  $(p_i - g_i)^2$'
    summary = (f'sum of squared errors = {s:.4f}\n'
               f'MSE = {s:.4f} / 4  =  {mse:.4f}')
    plot_contrib_bars(axes[1, i], x, sq, [str(t) for t in x],
                      formula, 'squared error', summary)
finalize(fig,
         'MSE breakdown — square each pointwise error, sum, divide by n',
         'method_mse_breakdown')


# =====================================================================
# 4. FRECHET  -- per-timepoint absolute deviation, take the MAX
#    (At n=4, optimal Frechet path is the diagonal for these genes,
#     so per-timepoint = per-path-step.)
# =====================================================================
def frechet_breakdown(p, g):
    devs = np.abs(p - g)
    return devs, float(devs.max())


genes_frechet = [('Gene A (match)', gene_A),
                 ('Gene E (small errs)', gene_E),
                 ('Gene F (1 disaster)', gene_F)]
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
for i, (name, g) in enumerate(genes_frechet):
    devs, fre = frechet_breakdown(pattern, g)
    max_idx = int(np.argmax(devs))
    plot_curves(axes[0, i], pattern, g, name, name)
    formula = r'contribution at $t_i$:  $|p_i - g_i|$'
    summary = f'Frechet = MAX of these  =  {fre:.4f}\n(red bar = winner)'
    plot_contrib_bars(axes[1, i], x, devs, [str(t) for t in x],
                      formula, '|deviation|', summary, max_idx=max_idx)
finalize(fig,
         'Frechet breakdown — take the largest single point-deviation along the optimal path',
         'method_frechet_breakdown')


# =====================================================================
# 3. DTW  -- per-timepoint sq-error along the OPTIMAL warping path
#    For genes where best path = diagonal, this matches MSE structure.
#    For Gene D (delayed spike), the optimal path warps non-trivially.
# =====================================================================
def dtw_breakdown(p, g):
    """Return (path_steps, step_costs, total) for the minimum-cost DTW path."""
    n, m = len(p), len(g)
    INF = float('inf')
    cost = np.full((n, m), INF)
    parent = {}
    cost[0, 0] = (p[0] - g[0]) ** 2
    for i in range(n):
        for j in range(m):
            if i == 0 and j == 0:
                continue
            d = (p[i] - g[j]) ** 2
            best, best_par = INF, None
            for di, dj in [(-1, -1), (-1, 0), (0, -1)]:
                if 0 <= i + di < n and 0 <= j + dj < m:
                    if cost[i + di, j + dj] + d < best:
                        best = cost[i + di, j + dj] + d
                        best_par = (i + di, j + dj)
            cost[i, j] = best
            parent[(i, j)] = best_par
    # Reconstruct path
    path = [(n - 1, m - 1)]
    while path[-1] != (0, 0):
        path.append(parent[path[-1]])
    path = path[::-1]
    step_costs = [(p[i] - g[j]) ** 2 for (i, j) in path]
    return path, step_costs, float(cost[n - 1, m - 1])


genes_dtw = [('Gene A (match)', gene_A),
             ('Gene C (wobble)', gene_C),
             ('Gene D (delayed spike)', gene_D)]
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
for i, (name, g) in enumerate(genes_dtw):
    path, step_costs, total = dtw_breakdown(pattern, g)
    # Top: curves with warping arrows
    ax = axes[0, i]
    ax.plot(x, pattern, 'o-', color='crimson', lw=2.2, label='Pattern', ms=7)
    ax.plot(x, g, 's-', color='steelblue', lw=2.2, label=name, ms=7)
    for (pi, gj) in path:
        ax.plot([x[pi], x[gj]], [pattern[pi], g[gj]],
                color='seagreen', alpha=0.6, lw=1.4)
    ax.set_xticks(x)
    ax.set_xlabel('Timepoint (h)')
    ax.set_ylim(-0.1, 1.2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title(f'{name}\n(green = DTW alignment path)', fontsize=10, fontweight='bold')

    # Bottom: per-step cost along the warping path
    ax = axes[1, i]
    step_pos = np.arange(len(step_costs))
    step_labels = [f'p{pi}\u2192g{gj}' for (pi, gj) in path]
    ax.bar(step_pos, step_costs, color='steelblue', edgecolor='black', linewidth=0.7)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xticks(step_pos)
    ax.set_xticklabels(step_labels, fontsize=8)
    ax.set_xlabel('warping path step')
    ax.set_ylabel('squared error')
    ax.set_title(r'per-step cost along path:  $(p_i - g_j)^2$', fontsize=9)
    ax.grid(alpha=0.3, axis='y')
    for xi, c in zip(step_pos, step_costs):
        ax.text(xi, c, f'{c:.4f}', ha='center', va='bottom', fontsize=7.5,
                fontweight='bold')
    summary = (f'sum along path = {total:.4f}\n'
               f'(DTW picks the path that minimizes this sum)')
    ax.text(0.5, -0.40, summary, transform=ax.transAxes,
            fontsize=9.5, ha='center', va='top', color='darkred', fontweight='bold',
            bbox=dict(facecolor='lightyellow', edgecolor='gray', boxstyle='round,pad=0.4'))

finalize(fig,
         'DTW breakdown — sum of squared errors along the optimal (warped) alignment path',
         'method_dtw_breakdown')


# =====================================================================
# 6. sMAPE  -- per-timepoint symmetric percentage error
# =====================================================================
def smape_breakdown(p, g):
    terms = []
    for pi, gi in zip(p, g):
        denom = abs(pi) + abs(gi)
        if denom == 0:
            terms.append(0.0)
        else:
            terms.append(2.0 * abs(pi - gi) / denom)
    terms = np.array(terms)
    return terms, float(terms.sum()), float(terms.mean())


# clusterThree pattern: zero at t=3 -> sMAPE failure case
p3 = np.array([1.0, 0.0, 0.91, 0.87])
gene_X = np.array([1.0, 0.01, 0.91, 0.87])  # near-perfect
gene_Y = np.array([0.5, 0.40, 0.60, 0.50])  # different shape
gene_Z = np.array([0.2, 0.90, 0.30, 0.40])  # OPPOSITE shape (peak at t=3)

genes_smape = [('Gene X (near-perfect)', gene_X),
               ('Gene Y (different)', gene_Y),
               ('Gene Z (opposite!)', gene_Z)]
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
for i, (name, g) in enumerate(genes_smape):
    terms, total, sm = smape_breakdown(p3, g)
    # Top: curves vs the clusterThree pattern
    ax = axes[0, i]
    ax.plot(x, p3, 'o-', color='crimson', lw=2.2,
            label='Pattern (zero at t=3)', ms=7)
    ax.plot(x, g, 's-', color='steelblue', lw=2.2, label=name, ms=7)
    ax.axvline(3, color='red', ls=':', alpha=0.4)
    ax.set_xticks(x)
    ax.set_xlabel('Timepoint (h)')
    ax.set_ylim(-0.1, 1.2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='center right')
    ax.set_title(name, fontsize=10, fontweight='bold')

    # Bottom: per-timepoint sMAPE term
    ax = axes[1, i]
    colors = ['red' if t >= 1.99 else 'steelblue' for t in terms]
    ax.bar(x, terms, color=colors, edgecolor='black', linewidth=0.7)
    ax.axhline(2.0, color='red', ls=':', alpha=0.6, lw=1.5)
    ax.set_xticks(x)
    ax.set_xlabel('Timepoint (h)')
    ax.set_ylabel('sMAPE term')
    ax.set_ylim(0, 2.35)
    ax.set_title(r'contribution at $t_i$:  $2|p_i - g_i| / (|p_i| + |g_i|)$',
                 fontsize=9)
    ax.grid(alpha=0.3, axis='y')
    for xi, c in zip(x, terms):
        ax.text(xi, c, f'{c:.3f}', ha='center', va='bottom', fontsize=8,
                fontweight='bold')
    summary = (f'sum = {total:.4f}\n'
               f'sMAPE = {total:.4f} / 4  =  {sm:.4f}\n'
               f'(red bar = saturated at 2.0)')
    ax.text(0.5, -0.40, summary, transform=ax.transAxes,
            fontsize=9.5, ha='center', va='top', color='darkred', fontweight='bold',
            bbox=dict(facecolor='lightyellow', edgecolor='gray', boxstyle='round,pad=0.4'))

finalize(fig,
         'sMAPE breakdown — average of percentage errors. Pattern zero at t=3 saturates every gene at 2.0',
         'method_smape_breakdown')

print('Wrote 6 breakdown figures to plots/method_*_breakdown.png')
