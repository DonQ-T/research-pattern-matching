"""Generate a 2x2 visual of the 4 constraint patterns for NOTES.md."""
import matplotlib.pyplot as plt

x = [0, 3, 6, 9]
patterns = {
    "clusterOne":   ([1.45,    12.0,   5.21,   6.51],   "Spike at t=3, partial recovery"),
    "clusterTwo":   ([0.371,   0.803,  1.22,   1.78],   "Steady monotonic increase"),
    "clusterThree": ([0.0405,  0.0,    0.0369, 0.0353], "Near-flat with dip at t=3"),
    "clusterFour":  ([0.00913, 0.638,  0.0178, 0.0243], "70x spike at t=3 (delta)"),
}

fig, axes = plt.subplots(2, 2, figsize=(10, 7))
for ax, (name, (vals, shape)) in zip(axes.flat, patterns.items()):
    ax.plot(x, vals, color="crimson", linewidth=2.5, marker="o", markersize=7)
    ax.set_title(f"{name} — {shape}", fontsize=11)
    ax.set_xlabel("Timepoint (h)")
    ax.set_ylabel("Expression")
    ax.set_xticks(x)
    ax.grid(alpha=0.3)
    for xi, yi in zip(x, vals):
        ax.annotate(f"{yi:g}", (xi, yi), textcoords="offset points",
                    xytext=(6, 6), fontsize=9, color="dimgray")

fig.suptitle("The 4 Constraint Patterns (raw values)", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig("plots/constraint_patterns_overview.png", dpi=130, bbox_inches="tight")
print("Wrote plots/constraint_patterns_overview.png")
