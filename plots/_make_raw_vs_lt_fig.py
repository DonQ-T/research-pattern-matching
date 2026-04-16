"""Side-by-side visual: raw constraint pattern vs log-transformed pattern,
one row per cluster."""
import matplotlib.pyplot as plt

x = [0, 3, 6, 9]

raw = {
    "clusterOne":   [1.45,    12.0,    5.21,    6.51],
    "clusterTwo":   [0.371,   0.803,   1.22,    1.78],
    "clusterThree": [0.0405,  0.0,     0.0369,  0.0353],
    "clusterFour":  [0.00913, 0.638,   0.0178,  0.0243],
}
lt = {
    "clusterOne":   [0.8945,  2.5687,  1.8260,  2.0165],
    "clusterTwo":   [0.0364,  0.5896,  0.7983,  1.0235],
    "clusterThree": [0.0397,  0.0,     0.0363,  0.0347],
    "clusterFour":  [0.00909, 0.4937,  0.01760, 0.02401],
}

fig, axes = plt.subplots(4, 2, figsize=(11, 11))

for i, name in enumerate(raw):
    # Raw column
    ax = axes[i, 0]
    vals = raw[name]
    ax.plot(x, vals, 'o-', color='crimson', lw=2.5, ms=8)
    ax.set_title(f'{name}  (raw)', fontsize=11, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xlabel('Timepoint (h)')
    ax.set_ylabel('Expression')
    ax.grid(alpha=0.3)
    for xi, yi in zip(x, vals):
        ax.annotate(f'{yi:g}', (xi, yi), textcoords='offset points',
                    xytext=(6, 6), fontsize=9, color='dimgray')
    nonzero = [v for v in vals if v > 0]
    if min(vals) == 0:
        ratio_label = 'has true zero'
    else:
        ratio_raw = max(vals) / min(nonzero)
        ratio_label = f'max/min = {ratio_raw:.0f}x'
    ax.text(0.97, 0.05, ratio_label,
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8.5, color='dimgray',
            bbox=dict(facecolor='white', edgecolor='lightgray', boxstyle='round,pad=0.3'))

    # LT column
    ax = axes[i, 1]
    vals = lt[name]
    ax.plot(x, vals, 'o-', color='steelblue', lw=2.5, ms=8)
    ax.set_title(f'{name}LT  (log-transformed)', fontsize=11, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xlabel('Timepoint (h)')
    ax.set_ylabel('Expression')
    ax.grid(alpha=0.3)
    for xi, yi in zip(x, vals):
        ax.annotate(f'{yi:g}', (xi, yi), textcoords='offset points',
                    xytext=(6, 6), fontsize=9, color='dimgray')
    nonzero = [v for v in vals if v > 0]
    if min(vals) == 0:
        ratio_label = 'has true zero'
    else:
        ratio_lt = max(vals) / min(nonzero)
        ratio_label = f'max/min = {ratio_lt:.0f}x'
    ax.text(0.97, 0.05, ratio_label,
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8.5, color='dimgray',
            bbox=dict(facecolor='white', edgecolor='lightgray', boxstyle='round,pad=0.3'))

fig.suptitle('Raw vs Log-Transformed Constraint Patterns\n'
             '(same shapes, but LT compresses the dynamic range)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('plots/constraint_patterns_raw_vs_lt.png', dpi=130, bbox_inches='tight')
print('Wrote plots/constraint_patterns_raw_vs_lt.png')
