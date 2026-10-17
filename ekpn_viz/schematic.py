"""Schematic of a single charged pore between two reservoirs (matplotlib, no data needed).

    python -m ekpn_viz.schematic docs/figures/single_pore_schematic.png
"""
import sys

import numpy as np


def single_pore_schematic(path="docs/figures/single_pore_schematic.png", dpi=200):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Ellipse, Polygon, Rectangle

    dark, light = "#245a8f", "#a9c8ea"
    cm = LinearSegmentedColormap.from_list("res", [dark, light])
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    ax.set_xlim(-0.2, 10.2); ax.set_ylim(-1.9, 2.0); ax.set_aspect("equal"); ax.axis("off")

    def cylinder(x0, x1, r, color, face):
        """Side-on cylinder: body + far cap (left) + visible face (right or left)."""
        ax.add_patch(Rectangle((x0, -r), x1 - x0, 2 * r, color=color, zorder=1))
        ax.add_patch(Ellipse((x0, 0), 0.7, 2 * r, color=color, zorder=1))
        ax.add_patch(Ellipse((x1, 0), 0.7, 2 * r, color=color, zorder=1))
        fx = x0 if face == "left" else x1
        ax.add_patch(Ellipse((fx, 0), 0.7, 2 * r, color="white", alpha=0.25, zorder=2))

    # reservoirs
    cylinder(0.3, 2.2, 1.5, dark, "left")
    cylinder(7.8, 9.7, 1.5, light, "right")
    # pore: horizontal gradient from dark to light
    h, x0, x1 = 0.55, 2.2, 7.8
    grad = np.linspace(0, 1, 256)[None, :]
    ax.imshow(grad, extent=(x0, x1, -h, h), cmap=cm, aspect="auto", zorder=3)
    for y in (h, -h):
        ax.plot([x0, x1], [y, y], "k--", lw=1, zorder=4)
    ax.add_patch(Ellipse((x1, 0), 0.25, 2 * h, fill=False, ec="k", lw=1.2, zorder=4))
    # EDL edges: thicker near the dilute side
    xx = np.linspace(x0, x1, 200)
    d = 0.12 + 0.28 * (xx - x0) / (x1 - x0)
    for s in (1, -1):
        ax.plot(xx, s * (h - d), "k-", lw=0.8, zorder=4)
    ax.annotate("", xy=(3.4, h), xytext=(3.4, h - 0.12 - 0.28 * 1.2 / 5.6),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, shrinkA=0, shrinkB=0,
                                mutation_scale=6), zorder=5)
    ax.text(3.4, 0.02, r"$\lambda_D$", ha="center", fontsize=11, zorder=5)
    # differential element
    for xe in (5.0, 5.35):
        ax.add_patch(Ellipse((xe, 0), 0.32, 2 * h + 0.25, fill=False, ec="#d62728", ls="--", lw=1.6, zorder=5))
    ax.annotate("", xy=(4.95, -h - 0.35), xytext=(5.4, -h - 0.35), arrowprops=dict(arrowstyle="<->", lw=0.9))
    ax.text(5.17, -h - 0.75, "dx", ha="center", fontsize=11)
    ax.text(5.17, h + 0.3, r"$\sigma$", ha="center", fontsize=14)
    # labels
    ax.text(1.3, 0.0, r"$C_{0L}$", color="white", ha="center", va="center", fontsize=15, zorder=3)
    ax.text(8.85, 0.25, r"$C_{0R}$", color="k", ha="center", va="center", fontsize=15, zorder=3)
    ax.text(8.85, -0.35, r"($C_{0L} > C_{0R}$)", ha="center", va="center", fontsize=9, zorder=3)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


if __name__ == "__main__":
    single_pore_schematic(*sys.argv[1:2])
