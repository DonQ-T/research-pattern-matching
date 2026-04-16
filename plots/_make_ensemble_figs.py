"""Visualize how the ensemble computation works:
   - Step 1: each gene has 6 ranks (one per method); votes and mean rank
             are computed from those ranks.
   - Step 2: the same set of (votes, mean) values can be sorted two ways,
             giving two different ensemble orderings.
"""
import numpy as np
import matplotlib.pyplot as plt

# -------- Toy example: 8 genes engineered to show the divergence --------
genes = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8']
descriptions = [
    'unanimous top',
    'loved by 4 methods, rejected by 2',
    'broadly agreed but mediocre',
    'unanimous solid 2nd-tier',
    'unanimous mid-tier',
    'top of 5 methods, 1 outlier',
    'bottom-of-top-20 across all',
    'never in any top-20',
]
methods = ['Pearson', 'DTW', 'Cosine', 'Frechet', 'MSE', 'sMAPE']
ranks = np.array([
    [ 2,  1,  3,  2,  1,  3],   # G1
    [ 1,  2,  1,  1, 35, 40],   # G2
    [18, 15, 19, 20, 16, 17],   # G3
    [ 4,  4,  4,  5,  4,  4],   # G4
    [ 8,  9,  7, 10,  8,  9],   # G5
    [50,  3,  2,  4,  5,  2],   # G6
    [12, 11, 18, 13, 19, 20],   # G7
    [100,105, 90,110, 85, 95],  # G8
])
votes = (ranks <= 20).sum(axis=1)
mean_ranks = ranks.mean(axis=1)

n_genes = len(genes)
n_methods = len(methods)

# ============================================================
# Figure 1 — the input: rank matrix + derived columns
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

display = np.clip(ranks, 1, 25).astype(float)
ax.imshow(display, cmap='RdYlGn_r', vmin=1, vmax=25, aspect='auto')

# Cell annotations
for i in range(n_genes):
    for j in range(n_methods):
        val = ranks[i, j]
        text_color = 'white' if display[i, j] > 16 else 'black'
        ax.text(j, i, str(val), ha='center', va='center',
                color=text_color, fontsize=10, fontweight='bold')

# Separator and side columns (votes, mean rank)
ax.axvline(n_methods - 0.5, color='black', linewidth=2)
for i in range(n_genes):
    if votes[i] == 6:
        col = '#4caf50'
    elif votes[i] >= 4:
        col = '#ffc107'
    else:
        col = '#f44336'
    ax.add_patch(plt.Rectangle((n_methods - 0.5, i - 0.5), 1, 1,
                               facecolor=col, alpha=0.35, edgecolor='black'))
    ax.text(n_methods, i, f'{votes[i]}', ha='center', va='center',
            fontsize=11, fontweight='bold')
    ax.add_patch(plt.Rectangle((n_methods + 0.5, i - 0.5), 1, 1,
                               facecolor='#9fa8da', alpha=0.35, edgecolor='black'))
    ax.text(n_methods + 1, i, f'{mean_ranks[i]:.1f}', ha='center', va='center',
            fontsize=11, fontweight='bold')

xticks = list(range(n_methods)) + [n_methods, n_methods + 1]
xlabels = methods + ['VOTES\n(/6)', 'MEAN\nRANK']
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, fontsize=10)
ax.set_yticks(range(n_genes))
ax.set_yticklabels([f'{g} — {d}' for g, d in zip(genes, descriptions)], fontsize=10)
ax.set_xlim(-0.5, n_methods + 1.5)
ax.set_title('Step 1 — each gene has 6 ranks (one per method).\n'
             'VOTES = count of ranks ≤ 20.    MEAN RANK = average of the 6 ranks.',
             fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('plots/ensemble_step1_ranks.png', dpi=130, bbox_inches='tight')
plt.close()


# ============================================================
# Figure 2 — same data, two different orderings
# ============================================================
order_votes = sorted(range(n_genes), key=lambda i: (-votes[i], mean_ranks[i]))
order_mean = sorted(range(n_genes), key=lambda i: mean_ranks[i])

fig, ax = plt.subplots(figsize=(12, 7))
ax.set_xlim(0, 14)
ax.set_ylim(-1.2, n_genes - 0.3)
ax.invert_yaxis()
ax.axis('off')

ax.text(3, -1, 'METHOD A — sort by VOTES\n(mean rank as tiebreak)',
        fontsize=11, fontweight='bold', ha='center', color='steelblue')
ax.text(11, -1, 'METHOD B — sort by MEAN RANK only',
        fontsize=11, fontweight='bold', ha='center', color='darkorange')

# Method A column
for pos, idx in enumerate(order_votes):
    label = f'{pos+1}.  {genes[idx]}   (votes={votes[idx]},  mean={mean_ranks[idx]:.1f})'
    ax.text(3, pos, label, ha='center', va='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#cce5ff',
                      edgecolor='steelblue', linewidth=1.5))

# Method B column
for pos, idx in enumerate(order_mean):
    label = f'{pos+1}.  {genes[idx]}   (votes={votes[idx]},  mean={mean_ranks[idx]:.1f})'
    ax.text(11, pos, label, ha='center', va='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffe5cc',
                      edgecolor='darkorange', linewidth=1.5))

# Connecting lines: red if position changes, gray otherwise
for idx in range(n_genes):
    pos_v = order_votes.index(idx)
    pos_m = order_mean.index(idx)
    if pos_v == pos_m:
        ax.plot([5.5, 8.5], [pos_v, pos_m], color='gray', alpha=0.35, linewidth=1)
    else:
        ax.plot([5.5, 8.5], [pos_v, pos_m], color='red', alpha=0.75, linewidth=2)

# Annotation for the genes that swap
ax.text(7, n_genes - 0.5,
        'G6 and G2 move UP under Method B\n'
        '(their good ranks pull the mean down\n'
        'even though some methods rejected them)\n\n'
        'G7 and G3 move DOWN under Method B\n'
        '(consistent mediocrity = high mean rank)',
        ha='center', va='top', fontsize=8.5,
        bbox=dict(facecolor='lightyellow', edgecolor='gray',
                  boxstyle='round,pad=0.5', alpha=0.9))

ax.set_title('Step 2 — same numbers, two ways to sort them.\n'
             'Red lines mark genes whose position changes between methods.',
             fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('plots/ensemble_step2_sorted.png', dpi=130, bbox_inches='tight')
plt.close()

print('Wrote plots/ensemble_step1_ranks.png and plots/ensemble_step2_sorted.png')
