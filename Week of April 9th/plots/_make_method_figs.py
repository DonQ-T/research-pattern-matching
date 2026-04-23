"""Generate per-method visual aids for NOTES.md."""
import numpy as np
import matplotlib.pyplot as plt

x = [0, 3, 6, 9]
xa = np.array(x)

# Reusable normalized vectors
pattern = np.array([0.0, 1.0, 0.356, 0.479])

gene_A = np.array([0.0, 1.0, 0.371, 0.486])  # shape match
gene_B = np.array([1.0, 1.0, 1.0, 1.0])      # housekeeping
gene_C = np.array([0.0, 1.0, 0.5, 0.8])      # wobble (over-stretched)
gene_D = np.array([0.0, 0.3, 1.0, 0.4])      # delayed spike at t=6
gene_E = np.array([0.05, 0.90, 0.40, 0.55])  # slightly off everywhere
gene_F = np.array([0.0, 1.0, 0.356, 0.9])    # perfect 3, bad at t=9


# ===============================================================
# 1. PEARSON — mean centering reveals shape
# ===============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
ax.plot(x, pattern, 'o-', color='crimson', linewidth=2.5, label='Pattern', markersize=8)
ax.plot(x, gene_A, 's--', color='steelblue', linewidth=2,
        label='Gene A (real shape match)', markersize=7)
ax.plot(x, gene_B, '^--', color='darkorange', linewidth=2,
        label='Gene B (flat housekeeping)', markersize=7)
ax.set_title('Raw normalized vectors', fontsize=11)
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Normalized expression')
ax.set_xticks(x)
ax.set_ylim(-0.6, 1.15)
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc='lower right')

ax = axes[1]
p_c = pattern - pattern.mean()
A_c = gene_A - gene_A.mean()
B_c = gene_B - gene_B.mean()
ax.axhline(0, color='black', linewidth=0.7)
ax.plot(x, p_c, 'o-', color='crimson', linewidth=2.5, label='Pattern (centered)', markersize=8)
ax.plot(x, A_c, 's--', color='steelblue', linewidth=2,
        label='Gene A:  r = 0.999', markersize=7)
ax.plot(x, B_c, '^--', color='darkorange', linewidth=2,
        label='Gene B:  collapses to zero (r = undef)', markersize=7)
ax.set_title('After mean-centering — what Pearson actually compares', fontsize=11)
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Centered expression')
ax.set_xticks(x)
ax.set_ylim(-0.6, 1.15)
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc='lower right')

fig.suptitle('Pearson: subtract the mean first → pure shape comparison',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/method_pearson.png', dpi=130, bbox_inches='tight')
plt.close()


# ===============================================================
# 2. COSINE — magnitude bias
# ===============================================================
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(x, pattern, 'o-', color='crimson', linewidth=2.5,
        label='Pattern (clusterOne, normalized)', markersize=9)
ax.plot(x, gene_A, 's--', color='steelblue', linewidth=2,
        label='Gene A (real match):   cos = 0.998   |  Pearson r = 0.999', markersize=7)
ax.plot(x, gene_C, 'D--', color='seagreen', linewidth=2,
        label='Gene C (wobble):         cos = 0.973   |  Pearson r = 0.95', markersize=7)
ax.plot(x, gene_B, '^--', color='darkorange', linewidth=2,
        label='Gene B (housekeeping):   cos = 0.788   |  Pearson r = undef', markersize=7)
ax.set_title('Cosine gives the flat housekeeping gene 0.79 — '
             'high enough to displace real matches',
             fontsize=11, fontweight='bold')
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Normalized expression')
ax.set_xticks(x)
ax.set_ylim(-0.05, 1.18)
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc='lower center')
ax.text(6.7, 0.27,
        "Without mean-centering, the\n"
        "constant vector [1,1,1,1] still\n"
        "has a positive angle with the\n"
        "pattern in 4-D space.",
        fontsize=8.5, color='dimgray', ha='center',
        bbox=dict(facecolor='lightyellow', edgecolor='gray', alpha=0.85))
plt.tight_layout()
plt.savefig('plots/method_cosine.png', dpi=130, bbox_inches='tight')
plt.close()


# ===============================================================
# 3. DTW — warping aligns time-shifted features
# ===============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Left: forced point-to-point (Euclidean / MSE) alignment
ax = axes[0]
ax.plot(x, pattern, 'o-', color='crimson', linewidth=2.5, label='Pattern (spike at t=3)', markersize=8)
ax.plot(x, gene_D, 's-', color='steelblue', linewidth=2.5,
        label='Gene D (spike at t=6)', markersize=8)
for xi, p, g in zip(x, pattern, gene_D):
    ax.plot([xi, xi], [p, g], color='gray', alpha=0.6, linewidth=1.5)
ax.set_title('Forced t-by-t alignment (Euclidean / MSE)\n'
             'MSE = 0.227  — peaks never line up', fontsize=10)
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Normalized expression')
ax.set_xticks(x)
ax.set_ylim(-0.05, 1.15)
ax.grid(alpha=0.3)
ax.legend(fontsize=8)

# Right: DTW warping path
ax = axes[1]
ax.plot(x, pattern, 'o-', color='crimson', linewidth=2.5, label='Pattern', markersize=8)
ax.plot(x, gene_D, 's-', color='steelblue', linewidth=2.5,
        label='Gene D', markersize=8)
# Optimal warping path (i,j): (0,0)→(1,1)→(1,2)→(2,3)→(3,3)
path = [(0, 0), (1, 1), (1, 2), (2, 3), (3, 3)]
for i, j in path:
    ax.plot([x[i], x[j]], [pattern[i], gene_D[j]],
            color='seagreen', alpha=0.75, linewidth=2)
ax.set_title('DTW: warp pattern\'s spike to align with gene\'s spike\n'
             'DTW distance ≈ 0.098  — 2.3× better than MSE', fontsize=10)
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Normalized expression')
ax.set_xticks(x)
ax.set_ylim(-0.05, 1.15)
ax.grid(alpha=0.3)
ax.legend(fontsize=8)

fig.suptitle('DTW: stretch the time axis to align matching features',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/method_dtw.png', dpi=130, bbox_inches='tight')
plt.close()


# ===============================================================
# 4. FRECHET — only the worst single point matters
# ===============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Gene E: small errors everywhere
ax = axes[0]
ax.plot(x, pattern, 'o-', color='crimson', linewidth=2.5, label='Pattern', markersize=8)
ax.plot(x, gene_E, 's-', color='steelblue', linewidth=2.5,
        label='Gene E (small errors everywhere)', markersize=8)
errs_E = (pattern - gene_E) ** 2
worst_E = int(np.argmax(errs_E))
for k, (xi, p, g) in enumerate(zip(x, pattern, gene_E)):
    color = 'red' if k == worst_E else 'gray'
    lw = 3 if k == worst_E else 2
    ax.plot([xi, xi], [p, g], color=color, linewidth=lw, alpha=0.85)
ax.set_title('Gene E: many small errors\n'
             'DTW=0.020   Frechet=0.010   (red = max)', fontsize=10)
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Normalized expression')
ax.set_xticks(x)
ax.set_ylim(-0.05, 1.15)
ax.grid(alpha=0.3)
ax.legend(fontsize=8)

# Gene F: perfect at 3, terrible at 1
ax = axes[1]
ax.plot(x, pattern, 'o-', color='crimson', linewidth=2.5, label='Pattern', markersize=8)
ax.plot(x, gene_F, 's-', color='steelblue', linewidth=2.5,
        label='Gene F (perfect 3, awful 1)', markersize=8)
errs_F = (pattern - gene_F) ** 2
worst_F = int(np.argmax(errs_F))
for k, (xi, p, g) in enumerate(zip(x, pattern, gene_F)):
    color = 'red' if k == worst_F else 'gray'
    lw = 3 if k == worst_F else 2
    ax.plot([xi, xi], [p, g], color=color, linewidth=lw, alpha=0.85)
ax.set_title('Gene F: one terrible point dominates\n'
             'DTW=0.177   Frechet=0.177   (worst pt = whole score)', fontsize=10)
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Normalized expression')
ax.set_xticks(x)
ax.set_ylim(-0.05, 1.15)
ax.grid(alpha=0.3)
ax.legend(fontsize=8)

fig.suptitle('Frechet: max error along the path  (vs DTW which sums)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/method_frechet.png', dpi=130, bbox_inches='tight')
plt.close()


# ===============================================================
# 5. MSE — average squared point-to-point error
# ===============================================================
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
genes = [
    ('Gene A (shape match)', gene_A, 'steelblue'),
    ('Gene C (wobble)', gene_C, 'seagreen'),
    ('Gene B (housekeeping)', gene_B, 'darkorange'),
]
for ax, (name, g, color) in zip(axes, genes):
    ax.plot(x, pattern, 'o-', color='crimson', linewidth=2.5, label='Pattern', markersize=8)
    ax.plot(x, g, 's-', color=color, linewidth=2, label=name, markersize=7)
    for xi, p, gi in zip(x, pattern, g):
        ax.fill_between([xi - 0.18, xi + 0.18], min(p, gi), max(p, gi),
                        color='gray', alpha=0.45)
    mse = float(np.mean((pattern - g) ** 2))
    ax.set_title(f'{name}\nMSE = {mse:.4f}', fontsize=10)
    ax.set_xlabel('Timepoint (h)')
    if ax is axes[0]:
        ax.set_ylabel('Normalized expression')
    ax.set_xticks(x)
    ax.set_ylim(-0.05, 1.18)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc='lower right')

fig.suptitle('MSE: gray rectangles = squared error at each timepoint, then averaged',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/method_mse.png', dpi=130, bbox_inches='tight')
plt.close()


# ===============================================================
# 6. sMAPE — failure on a pattern with a zero
# ===============================================================
# clusterThree pattern (normalized): near-flat with a zero/near-zero at t=3
p3 = np.array([1.0, 0.0, 0.91, 0.87])

# Three genes with very different shapes — but all with a nonzero at t=3
gene_X = np.array([1.0, 0.01, 0.91, 0.87])   # nearly perfect, tiny noise at t=3
gene_Y = np.array([0.5, 0.40, 0.60, 0.50])   # very different shape
gene_Z = np.array([0.2, 0.90, 0.30, 0.40])   # *opposite* shape (peak at t=3!)


def smape_term(p, g):
    if p == 0 and g == 0:
        return 0.0
    return 2.0 * abs(p - g) / (abs(p) + abs(g))


fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

ax = axes[0]
ax.plot(x, p3, 'o-', color='crimson', linewidth=2.5,
        label='clusterThree pattern (zero at t=3)', markersize=9)
ax.plot(x, gene_X, 's--', color='seagreen', linewidth=2,
        label='Gene X (near-perfect)', markersize=7)
ax.plot(x, gene_Y, 'D--', color='steelblue', linewidth=2,
        label='Gene Y (different shape)', markersize=7)
ax.plot(x, gene_Z, '^--', color='darkorange', linewidth=2,
        label='Gene Z (opposite shape!)', markersize=7)
ax.axvline(3, color='red', linestyle=':', alpha=0.5)
ax.text(3.1, 0.45, 'pattern = 0 here', color='red', fontsize=9)
ax.set_title('Three very different genes vs clusterThree pattern', fontsize=10)
ax.set_xlabel('Timepoint (h)')
ax.set_ylabel('Normalized expression')
ax.set_xticks(x)
ax.set_ylim(-0.05, 1.15)
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc='center right')

# Per-timepoint sMAPE bar chart
ax = axes[1]
terms_X = [smape_term(p, g) for p, g in zip(p3, gene_X)]
terms_Y = [smape_term(p, g) for p, g in zip(p3, gene_Y)]
terms_Z = [smape_term(p, g) for p, g in zip(p3, gene_Z)]
xpos = np.arange(4)
w = 0.27
ax.bar(xpos - w, terms_X, w, color='seagreen',
       label=f'Gene X  (mean = {np.mean(terms_X):.3f})')
ax.bar(xpos,     terms_Y, w, color='steelblue',
       label=f'Gene Y  (mean = {np.mean(terms_Y):.3f})')
ax.bar(xpos + w, terms_Z, w, color='darkorange',
       label=f'Gene Z  (mean = {np.mean(terms_Z):.3f})')
ax.axhline(2.0, color='red', linestyle=':', alpha=0.6, linewidth=1.5,
           label='sMAPE saturation (2.0)')
ax.set_xticks(xpos)
ax.set_xticklabels(['t=0', 't=3', 't=6', 't=9'])
ax.set_ylabel('sMAPE contribution at this timepoint')
ax.set_title('At t=3 every gene scores exactly 2.0 → no signal\n'
             'Ranking is decided by the other 3 (uninformative) timepoints',
             fontsize=10)
ax.set_ylim(0, 2.35)
ax.grid(alpha=0.3, axis='y')
ax.legend(fontsize=8, loc='upper right')

fig.suptitle('sMAPE: a true zero in the pattern saturates the metric',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/method_smape.png', dpi=130, bbox_inches='tight')
plt.close()

print('Wrote 6 method figures to plots/method_*.png')
