"""Generic plots that work for any network: network map, path profiles, per-pore grid, currents."""
import os

import numpy as np

from .io import FIELD_LABELS, read_current_flow, read_log
from .layout import node_positions


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    return plt, animation


def save_animation(fig, update, n_frames, base, gif=True, mp4=False, fps=14, gif_dpi=70, mp4_dpi=130,
                   gif_max_frames=60):
    """Write `base`.gif (subsampled to <= gif_max_frames, low dpi, small file) and/or `base`.mp4
    (every frame, H.264 via ffmpeg). `update(k)` draws frame k. Returns the written paths."""
    _, animation = _mpl()
    paths = []
    if mp4:
        if not animation.writers.is_available("ffmpeg"):
            raise SystemExit("--mp4 needs ffmpeg on PATH (apt install ffmpeg, or pip install imageio-ffmpeg "
                             "and set matplotlib rcParams['animation.ffmpeg_path'])")
        w = animation.FFMpegWriter(fps=fps, codec="libx264", bitrate=-1,
                                   extra_args=["-pix_fmt", "yuv420p", "-crf", "23"])
        animation.FuncAnimation(fig, update, frames=n_frames).save(base + ".mp4", writer=w, dpi=mp4_dpi)
        paths.append(base + ".mp4")
    if gif:
        step = max(1, int(np.ceil(n_frames / gif_max_frames)))
        ks = list(range(0, n_frames, step))
        if ks[-1] != n_frames - 1:
            ks.append(n_frames - 1)
        animation.FuncAnimation(fig, update, frames=ks).save(
            base + ".gif", writer=animation.PillowWriter(fps=max(4, fps * len(ks) // max(n_frames, 1) or 8)),
            dpi=gif_dpi)
        paths.append(base + ".gif")
    return paths


def _slit_geometry(run, pos, width):
    """Flat rectangular channels (slit pores) for a 2D schematic.
    Returns per-cell quads (coloured by the field), per-pore outlines, and the half-widths.
    The drawn width is `width` (plot units) times a factor from the pore area (area / max area,
    floored at 0.5 so thin elements stay visible), mimicking slit channels of different height."""
    pores = run.net["pores"]
    areas = np.array([float(pr.get("area", 1.0)) for pr in pores])
    amax = areas.max() if len(areas) else 1.0
    half = [0.5 * width * max(0.5, a_ / amax) for a_ in areas]
    quads, outlines = [], []
    for p, pr in enumerate(pores):
        a, b = pos[pr["ends"][0]], pos[pr["ends"][1]]
        d = b - a
        L = np.linalg.norm(d) + 1e-12
        n = np.array([-d[1], d[0]]) / L
        hw = half[p]
        # a thinner pore starts at the edge of the wider channels it joins (drawn inset), so its
        # whole coloured length stays visible instead of disappearing under the channel slabs
        ins = []
        for e in pr["ends"]:
            wider = [half[q] for q, pq in enumerate(pores) if q != p and e in pq["ends"] and half[q] > hw]
            ins.append(min(max(wider, default=0.0) / L, 0.45))
        a, b = a + ins[0] * d, b - ins[1] * d
        d = b - a
        m = len(run.field(0, p, "C_bar"))
        s_ = np.linspace(0, 1, m + 1)
        for i in range(m):
            p0, p1 = a + s_[i] * d, a + s_[i + 1] * d
            quads.append([p0 - hw * n, p1 - hw * n, p1 + hw * n, p0 + hw * n])
        outlines.append([a - hw * n, b - hw * n, b + hw * n, a + hw * n])
    return quads, outlines, np.array(half)


def _oblique(xy, shear=0.55, squash=0.5):
    """Tilted-plane view: x' = x + shear*y, y' = squash*y (keeps straight lines straight)."""
    xy = np.asarray(xy, float)
    return np.stack([xy[..., 0] + shear * xy[..., 1], squash * xy[..., 1]], axis=-1)


def _draw_slits(ax, run, pos, k, field, norm, cmap, width, depth=0.45, oblique=True):
    """Draw slit channels as flat slabs seen from above at an angle (like the README schematic):
    coloured top faces, grey side walls of thickness `depth*width`. Returns the PolyCollection
    whose array holds the per-cell field (update with set_array)."""
    from matplotlib.collections import PolyCollection
    quads, outlines, half = _slit_geometry(run, pos, width)
    span = np.ptp(np.asarray(pos, float), axis=0)
    # fixed oblique view, like the README schematic; layouts with almost no vertical extent
    # (e.g. a single row of pores) get a little more squash so the slabs remain readable
    ratio = float(span[0] / max(span[1], 1e-9))
    shear, squash = 0.6, (0.6 if ratio < 8 else 1.0)
    T = (lambda a: _oblique(a, shear, squash)) if oblique else (lambda a: np.asarray(a, float))
    quads = [T(q) for q in quads]
    outlines = [T(o) for o in outlines]
    pos = T(pos)
    off = np.array([0.0, -1.0]) * width * depth
    shade = [np.array(o) + off for o in outlines]
    ax.add_collection(PolyCollection(shade, facecolors="0.55", edgecolors="0.35", linewidths=0.6, zorder=1))
    for o in outlines:   # side walls joining top face and shadow
        o = np.array(o)
        for i in range(4):
            q = [o[i], o[(i + 1) % 4], o[(i + 1) % 4] + off, o[i] + off]
            ax.add_collection(PolyCollection([q], facecolors="0.7", edgecolors="0.4", linewidths=0.4, zorder=1.5))
    pc = PolyCollection(quads, cmap=cmap, norm=norm, edgecolors="none", zorder=2)
    pc.set_array(np.concatenate([run.field(k, p, field) for p in range(run.n_pores)]))
    ax.add_collection(pc)
    ax.add_collection(PolyCollection(outlines, facecolors="none", edgecolors="0.15", linewidths=0.8, zorder=3))
    res = run.net.get("reservoirs") or []
    hw = half.max() if len(half) else width / 2
    for nd in range(run.n_nodes):
        fixed = bool(res) and nd < len(res) and res[nd]["potential_known"]
        c = pos[nd]
        sq = [c + np.array(v) for v in T(hw * np.array([(-1, -1), (1, -1), (1, 1), (-1, 1)], float))]
        ax.add_collection(PolyCollection([sq], facecolors="0.15" if fixed else "0.85",
                                         edgecolors="0.1", linewidths=0.8, zorder=4))
        if fixed:
            ax.annotate(f"V={res[nd]['potential']:g}", c, textcoords="offset points", xytext=(0, 10),
                        ha="center", fontsize=8, zorder=5)
    return pc


def _norm(run, field, values):
    """Colour scale: C_bar diverging around the initial bulk value (blue depleted, red
    enriched); signed fields (U, P, Mu) diverging around 0; else linear."""
    from matplotlib.colors import Normalize, TwoSlopeNorm
    lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
    if field in ("C_bar", "C"):
        bulk = float(np.median(np.concatenate([run.field(0, p, field) for p in range(run.n_pores)])))
        lo = min(lo, 0.0)
        if lo < bulk < hi:
            return TwoSlopeNorm(vmin=lo, vcenter=bulk, vmax=hi), "RdBu_r"
        return Normalize(lo, hi), "viridis"
    if field == "U" or lo < 0 < hi:   # signed: diverging, centred at 0
        m = max(-lo, hi) or 1.0
        return TwoSlopeNorm(vmin=-m, vcenter=0, vmax=m), "PuOr_r"
    return Normalize(lo, hi if hi > lo else lo + 1), "viridis"


def _label(field):
    return FIELD_LABELS.get(field, field)


def _pick_frames(run, times):
    if times is None:
        n = len(run.times)
        idx = sorted({0, n // 10, n // 3, n - 1})
    else:
        idx = [int(np.argmin(abs(run.times - t))) for t in times]
    return idx


# ---------------- network map ----------------
def pore_means(run, k, field="U"):
    """Mean of a field over each pore's cells at frame k (for U: net flow direction)."""
    return np.array([float(np.mean(run.field(k, p, field))) for p in range(run.n_pores)])


def simple_cycles(n_nodes, edges, max_len=12, max_count=2000):
    """All simple cycles of an undirected multigraph, as lists of (pore, from, to) steps."""
    adj = {i: [] for i in range(n_nodes)}
    for p, (a, b) in enumerate(edges):
        adj[a].append((p, b)); adj[b].append((p, a))
    found, seen = [], set()

    def dfs(start, node, path, used, visited):
        if len(found) >= max_count or len(path) > max_len:
            return
        for p, nb in adj[node]:
            if p in used:
                continue
            if nb == start and len(path) >= 1:
                key = frozenset(used | {p})
                if key not in seen and (len(path) >= 2 or edges[p] != edges[path[0][0]] or p != path[0][0]):
                    seen.add(key); found.append(path + [(p, node, nb)])
            elif nb not in visited and nb > start:
                dfs(start, nb, path + [(p, node, nb)], used | {p}, visited | {nb})

    for s in range(n_nodes):
        dfs(s, s, [], frozenset(), {s})
    return found


def circulating_loops(run, k, rel_tol=1e-3):
    """Cycles where the net flow (mean U per pore) runs the same way all around the loop.
    Returns list of dicts(pores, nodes, strength=min |mean U| along the loop)."""
    u = pore_means(run, k, "U")
    tol = rel_tol * (np.abs(u).max() or 1.0)
    edges = [tuple(p["ends"]) for p in run.net["pores"]]
    loops = []
    for cyc in simple_cycles(run.n_nodes, edges):
        signs = []
        for p, a, b in cyc:
            along = u[p] if edges[p][0] == a else -u[p]   # flow measured in the walking direction
            signs.append(along)
        signs = np.array(signs)
        if np.all(np.abs(signs) > tol) and (np.all(signs > 0) or np.all(signs < 0)):
            d = 1 if signs[0] > 0 else -1
            nodes = [cyc[0][1]] + [b for _, _, b in cyc]
            if d < 0:
                nodes = nodes[::-1]
            loops.append(dict(pores=[p for p, _, _ in cyc], nodes=nodes, strength=float(np.abs(signs).min())))
    return loops


def network_map(run, field="C_bar", times=None, out="docs/figures", prefix="network", layout=None,
                width=6, title=None, arrows=None, loops=None, shape="line", slit_width=None):
    """Schematic of the whole network, every pore drawn node-to-node and coloured cell by cell.
    Pore lengths are not to scale. Squares: fixed-potential reservoirs; circles: internal nodes."""
    plt, _ = _mpl()
    from matplotlib.collections import LineCollection
    pos, src = node_positions(run, layout)
    idx = _pick_frames(run, times)
    allv = np.concatenate([run.field(k, p, field) for k in idx for p in range(run.n_pores)])
    norm, cmap = _norm(run, field, allv)
    ncol = min(len(idx), 1 if shape == "slit" else 2)
    nrow = int(np.ceil(len(idx) / ncol))
    span = np.ptp(pos, axis=0)
    aspect = (span[1] + 0.3) / (span[0] + 0.3)
    if shape == "slit":
        fig, axs = plt.subplots(nrow, ncol, figsize=(8.5, 1.7 * nrow + 0.7), squeeze=False)
    else:
        fig, axs = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, max(2.2, 5.2 * aspect) * nrow + 0.6), squeeze=False)
    res = run.net.get("reservoirs") or []
    arrows = (field == "U") if arrows is None else arrows
    loops = (field == "U") if loops is None else loops
    slit_w = slit_width or 0.06 * float(np.ptp(pos, axis=0).max() or 1.0)
    for ax, k in zip(axs.flat, idx):
        if shape == "slit":
            lc = _draw_slits(ax, run, pos, k, field, norm, cmap, slit_w)
            ax.autoscale(); ax.margins(0.08); ax.set_aspect("equal"); ax.axis("off")
            ax.set_title(f"t = {run.times[k]:.3g}", fontsize=10)
            continue
        segs, vals = [], []
        for p, pore in enumerate(run.net["pores"]):
            a, b = pos[pore["ends"][0]], pos[pore["ends"][1]]
            v = run.field(k, p, field)
            s = np.linspace(0, 1, len(v) + 1)[:, None]
            pts = a + s * (b - a)
            segs += [[pts[i], pts[i + 1]] for i in range(len(v))]
            vals.append(v)
        lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=width, capstyle="butt")
        lc.set_array(np.concatenate(vals))
        ax.add_collection(lc)
        for nd in range(run.n_nodes):
            fixed = bool(res) and nd < len(res) and res[nd]["potential_known"]
            ax.plot(*pos[nd], "s" if fixed else "o", ms=8 if fixed else 6,
                    color="k" if fixed else "w", mec="k", zorder=3)
            lab = f"{nd}" + (f"\nV={res[nd]['potential']:g}" if fixed else "")
            ax.annotate(lab, pos[nd], textcoords="offset points", xytext=(6, 6), fontsize=7)
        for p, pore in enumerate(run.net["pores"]):
            m = (pos[pore["ends"][0]] + pos[pore["ends"][1]]) / 2
            ax.annotate(f"p{p}", m, textcoords="offset points", xytext=(0, -12), fontsize=7,
                        ha="center", color="0.35")
        ttl = f"t = {run.times[k]:.3g}"
        if arrows:   # net flow direction per pore, arrow size ~ |mean U|
            um = pore_means(run, k, "U")
            umax = np.abs(um).max() or 1.0
            for p, pore in enumerate(run.net["pores"]):
                a, b = pos[pore["ends"][0]], pos[pore["ends"][1]]
                if abs(um[p]) < 1e-3 * umax:
                    continue
                d = (b - a) / (np.linalg.norm(b - a) + 1e-12) * np.sign(um[p])
                c = (a + b) / 2
                L = 0.12 * np.ptp(pos, axis=0).max() * (0.4 + 0.6 * abs(um[p]) / umax)
                ax.annotate("", xy=c + d * L / 2, xytext=c - d * L / 2,
                            arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2, mutation_scale=12), zorder=4)
        if loops:
            found = circulating_loops(run, k)
            for lp in found:
                poly = np.array([pos[n] for n in lp["nodes"]])
                ax.fill(poly[:, 0], poly[:, 1], color="gold", alpha=0.25, zorder=0)
                cen = poly[:-1].mean(0)
                ax.text(*cen, "↻" if _orientation(poly) < 0 else "↺", fontsize=18, ha="center",
                        va="center", color="darkgoldenrod", zorder=5)
            if found:
                ttl += f"  ({len(found)} circulating loop{'s' * (len(found) > 1)})"
        ax.autoscale(); ax.margins(0.12); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(ttl, fontsize=10)
    for ax in list(axs.flat)[len(idx):]:
        ax.axis("off")
    fig.colorbar(lc, ax=axs, shrink=0.8, label=_label(field))
    fig.suptitle(title or f"{_label(field)} on the network (layout: {src}; pore lengths not to scale)",
                 fontsize=10)
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"{prefix}_{field}_{'slit' if shape == 'slit' else 'map'}.png")
    fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)
    return path


def _orientation(poly):
    """Signed area (>0: counter-clockwise)."""
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])


def flow_report(run, out="docs/figures", prefix="flow", rel_tol=1e-3):
    """Mean U per pore vs time with sign reversals marked, plus a text summary of
    reversals, in-pore stagnation points (U changes sign along a pore) and circulating loops."""
    plt, _ = _mpl()
    t = run.times
    U = np.array([pore_means(run, k, "U") for k in range(len(t))])
    tol = rel_tol * (np.abs(U).max() or 1.0)
    Us = np.where(np.abs(U) > tol, np.sign(U), 0)   # sign, 0 = "no net flow" (below tolerance)

    def flips_of(p):
        nz = np.nonzero(Us[:, p])[0]
        return [nz[i + 1] for i in range(len(nz) - 1) if Us[nz[i], p] != Us[nz[i + 1], p]]
    fig, ax = plt.subplots(figsize=(8, 4.4))
    lines = []
    for p in range(run.n_pores):
        a, b = run.net["pores"][p]["ends"]
        ln, = ax.plot(t, U[:, p], lw=1.6, label=f"p{p} ({a}→{b})")
        lines.append(ln)
        fl = flips_of(p)
        ax.plot(t[fl], U[fl, p], "o", ms=7, mfc="none", color=ln.get_color())
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel("t (nondimensional)"); ax.set_ylabel("mean U in pore (+ = first end → second)")
    ax.set_title("Net flow direction per pore (circles: sign reversals)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"{prefix}_direction.png")
    fig.tight_layout(); fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)

    rep = []
    for p in range(run.n_pores):
        fl = flips_of(p)
        if fl:
            rep.append(f"pore {p}: net flow reverses at t = " + ", ".join(f"{t[i]:.4g}" for i in fl))
        if not Us[:, p].any():
            rep.append(f"pore {p}: no net flow (|mean U| < {tol:.2g} throughout)")
    k = len(t) - 1
    for p in range(run.n_pores):
        u = run.field(k, p, "U")
        if np.abs(u).max() < tol:
            continue
        sc = np.nonzero(np.sign(u[1:]) * np.sign(u[:-1]) < 0)[0]
        if len(sc):
            rep.append(f"pore {p}: U changes sign inside the pore at the final time ({len(sc)} stagnation point(s))")
    first = None
    for kk in range(len(t)):
        if circulating_loops(run, kk):
            first = kk
            break
    final = circulating_loops(run, k)
    for lp in final:
        rep.append("circulating loop at final time: nodes " + " → ".join(map(str, lp["nodes"])) +
                   f" (pores {lp['pores']}, min |mean U| on loop {lp['strength']:.3g})")
    if first is not None:
        rep.append(f"circulation first detected at t = {t[first]:.4g}")
    if not final and first is None:
        rep.append("no circulating loop detected")
    return path, rep


# ---------------- path along chosen pores ----------------
def path_plots(run, pores, field="C_bar", threshold=None, out="docs/figures", prefix="path",
               gif=True, title=None, mp4=False, log_time=False):
    """Field along a user-chosen chain of pores: space-time map with tracked fronts, profiles,
    front positions vs time, and an optional GIF. A front is every place where the field
    crosses `threshold` (default: half the initial median for C_bar, else not tracked)."""
    plt, animation = _mpl()
    x, _, joints, seq = run.path(pores, field, 0)
    V = np.array([run.path(pores, field, k)[1] for k in range(len(run.times))])
    t = run.times
    if threshold is None and field in ("C_bar", "C"):
        threshold = 0.5 * float(np.median(V[0]))
    norm, cmap = _norm(run, field, V)
    os.makedirs(out, exist_ok=True)
    pre = os.path.join(out, f"{prefix}_{field}")
    name = title or ("path " + " → ".join(f"p{p}" for p in pores) + "  (nodes " + "-".join(map(str, seq)) + ")")
    tl = "t (nondimensional)"

    fronts = []   # (k, x) of every crossing
    if threshold is not None:
        for k in range(len(t)):
            s = np.sign(V[k] - threshold)
            ii = np.nonzero(s[1:] * s[:-1] < 0)[0]
            for i in ii:
                xf = x[i] + (threshold - V[k][i]) * (x[i + 1] - x[i]) / (V[k][i + 1] - V[k][i])
                fronts.append((k, xf))
    fr = np.array(fronts) if fronts else np.zeros((0, 2))

    def junctions(ax):
        for j, nd in zip(joints, seq[1:-1]):
            ax.axvline(j, color="k", lw=0.8, ls=":")
            ax.text(j, 1.0, f" node {nd}", transform=ax.get_xaxis_transform(), fontsize=7, va="bottom")

    fig, ax = plt.subplots(figsize=(8, 4.8))
    m = ax.pcolormesh(x, t, V, shading="auto", cmap=cmap, norm=norm, rasterized=True)
    if len(fr):
        ax.plot(fr[:, 1], t[fr[:, 0].astype(int)], "k.", ms=1.5, label=f"front ({_label(field)} = {threshold:.3g})")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.5)
    junctions(ax)
    if log_time and len(t) > 1:
        ax.set_yscale("log"); ax.set_ylim(t[t > 0].min(), t[-1])
    ax.set_xlabel("x along path"); ax.set_ylabel(tl + (" (log scale)" if log_time else ""))
    ax.set_title(name, fontsize=10, pad=14)
    fig.colorbar(m, ax=ax, label=_label(field))
    fig.tight_layout(); fig.savefig(pre + "_spacetime.png", dpi=160); plt.close(fig)

    sel = np.unique(np.r_[0, np.geomspace(1, len(t) - 1, 10).astype(int)]) if len(t) > 1 else [0]
    colors = plt.get_cmap("plasma")(np.linspace(0, 0.9, len(sel)))
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for c_, k in zip(colors, sel):
        ax.plot(x, V[k], color=c_, lw=1.5, label=f"t = {t[k]:.3g}")
    if threshold is not None:
        ax.axhline(threshold, color="0.5", lw=0.8, ls="--")
    junctions(ax)
    ax.set_xlabel("x along path"); ax.set_ylabel(_label(field)); ax.set_title(name, fontsize=10, pad=14)
    ax.legend(fontsize=8, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout(); fig.savefig(pre + "_profiles.png", dpi=160, bbox_inches="tight"); plt.close(fig)

    if len(fr):
        fig, ax = plt.subplots(figsize=(8, 4.2))
        ax.plot(t[fr[:, 0].astype(int)], fr[:, 1], ".", ms=2)
        for j in joints:
            ax.axhline(j, color="k", lw=0.8, ls=":")
        if log_time:
            ax.set_xscale("log")
        ax.set_xlabel(tl); ax.set_ylabel("front position x along path")
        ax.set_ylim(0, x[-1] + (x[1] - x[0]) / 2)
        ax.set_title(f"fronts where {_label(field)} = {threshold:.3g}", fontsize=10)
        fig.tight_layout(); fig.savefig(pre + "_fronts.png", dpi=160); plt.close(fig)

    if gif or mp4:
        fig, ax = plt.subplots(figsize=(6.6, 3.8))
        ln, = ax.plot(x, V[0], lw=2, color="#1f77b4")
        lo, hi = np.nanmin(V), np.nanmax(V)
        ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.08 * (hi - lo))
        junctions(ax)
        ax.set_xlabel("x along path"); ax.set_ylabel(_label(field))
        ttl = ax.set_title("")
        fig.tight_layout()

        def upd(k):
            ln.set_ydata(V[k]); ttl.set_text(f"{name}\nt = {t[k]:.4f}")
            return ln, ttl

        save_animation(fig, upd, len(t), pre, gif=gif, mp4=mp4)
        plt.close(fig)
    return dict(x=x, V=V, fronts=fr, nodes=seq)


# ---------------- per-pore grid ----------------
def pores_grid(run, field="C_bar", out="docs/figures", prefix="pores", pores=None):
    """Small-multiple space-time maps, one per pore, on a shared colour scale."""
    plt, _ = _mpl()
    pores = list(range(run.n_pores)) if pores is None else pores
    data = {p: np.array([run.field(k, p, field) for k in range(len(run.times))]) for p in pores}
    norm, cmap = _norm(run, field, np.concatenate([d.ravel() for d in data.values()]))
    ncol = min(4, len(pores))
    nrow = int(np.ceil(len(pores) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.8 * nrow), squeeze=False, layout="constrained")
    for ax, p in zip(axs.flat, pores):
        m = ax.pcolormesh(run.x(p), run.times, data[p], shading="auto", cmap=cmap, norm=norm, rasterized=True)
        a, b = run.net["pores"][p]["ends"]
        ax.set_title(f"pore {p}: node {a} → {b}", fontsize=9)
        ax.set_xlabel("x", fontsize=8); ax.set_ylabel("t", fontsize=8); ax.tick_params(labelsize=7)
    for ax in list(axs.flat)[len(pores):]:
        ax.axis("off")
    fig.colorbar(m, ax=axs, shrink=0.8, label=_label(field))
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"{prefix}_{field}_grid.png")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


# ---------------- currents ----------------
def currents(run, out="docs/figures", prefix="currents"):
    """I_tot(t) from the log (if given) and the final per-pore current I and flow Q
    (current_flow_area.bin). Sign: positive = from a pore's first end to its second."""
    plt, _ = _mpl()
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.8))
    if run.log:
        tl, it = read_log(run.log)
        k = tl <= run.times[-1] + 1e-12
        axs[0].plot(tl[k], it[k], "C3-", lw=1.5)
        axs[0].set_xlabel("t (nondimensional)"); axs[0].set_ylabel(r"$I_{\rm tot}$")
        axs[0].set_title("total current (solver log)", fontsize=10)
    else:
        axs[0].axis("off")
    cf_path = os.path.join(run.run_dir, "current_flow_area.bin")
    cf = None
    if os.path.exists(cf_path):
        cf = read_current_flow(cf_path)
        idx = np.arange(len(cf["I"]))
        axs[1].bar(idx - 0.2, cf["I"], 0.4, label="current I")
        axs[1].bar(idx + 0.2, cf["Q"], 0.4, label="flow rate Q")
        axs[1].axhline(0, color="k", lw=0.6)
        axs[1].set_xticks(idx, [f"p{i}\n{run.net['pores'][i]['ends'][0]}→{run.net['pores'][i]['ends'][1]}"
                                for i in idx], fontsize=8)
        axs[1].set_title("per-pore current and flow, final time (+ = first end → second)", fontsize=10)
        axs[1].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"{prefix}.png")
    fig.savefig(path, dpi=160); plt.close(fig)
    return path, cf


def node_balance(run, cf):
    """Kirchhoff check at every node: sum of currents in minus out (should be ~0 at internal nodes)."""
    bal = np.zeros(run.n_nodes)
    for p, pore in enumerate(run.net["pores"]):
        a, b = pore["ends"]
        bal[a] -= cf["I"][p]; bal[b] += cf["I"][p]
    return bal


# ---------------- network + paths animation ----------------
def animate(run, paths, field="C_bar", out="docs/figures", prefix="anim", layout=None, gif=True, mp4=True,
            fps=12, title=None, width=6, shape="line"):
    """Animation: network schematic coloured by `field` (left/top) plus the field along each chain
    of pores in `paths` (list of pore lists). Uses the snapshot times as frames (uneven spacing is
    kept, so early fast dynamics get more frames when the run wrote them)."""
    plt, _ = _mpl()
    from matplotlib.collections import LineCollection
    pos, src = node_positions(run, layout)
    nT = len(run.times)
    allv = np.concatenate([run.field(k, p, field) for k in range(nT) for p in range(run.n_pores)])
    norm, cmap = _norm(run, field, allv)
    span = np.ptp(pos, axis=0)
    fig = plt.figure(figsize=(11, 4.6))
    gs = fig.add_gridspec(len(paths), 2, width_ratios=[1.05, 1.25], wspace=0.25, hspace=0.55)
    axn = fig.add_subplot(gs[:, 0])
    res = run.net.get("reservoirs") or []
    if shape == "slit":
        lc = _draw_slits(axn, run, pos, 0, field, norm, cmap, 0.06 * float(np.ptp(pos, axis=0).max() or 1.0))
    else:
        segs = []
        for p, pore in enumerate(run.net["pores"]):
            a, b = pos[pore["ends"][0]], pos[pore["ends"][1]]
            n = len(run.field(0, p, field))
            s_ = np.linspace(0, 1, n + 1)[:, None]
            pts = a + s_ * (b - a)
            segs += [[pts[i], pts[i + 1]] for i in range(n)]
        lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=width, capstyle="butt")
        lc.set_array(np.concatenate([run.field(0, p, field) for p in range(run.n_pores)]))
        axn.add_collection(lc)
    for nd in (range(run.n_nodes) if shape != "slit" else []):
        fixed = bool(res) and nd < len(res) and res[nd]["potential_known"]
        axn.plot(*pos[nd], "s" if fixed else "o", ms=8 if fixed else 6, color="k" if fixed else "w", mec="k", zorder=3)
        if fixed:
            axn.annotate(f"V={res[nd]['potential']:g}", pos[nd], textcoords="offset points", xytext=(0, 8),
                         ha="center", fontsize=8)
    axn.autoscale(); axn.margins(0.08); axn.set_aspect("equal"); axn.axis("off")
    axn.set_anchor("C")
    fig.colorbar(lc, ax=axn, orientation="vertical", fraction=0.04, pad=0.02, shrink=0.6, label=_label(field))
    lines, data = [], []
    for i, pores in enumerate(paths):
        ax = fig.add_subplot(gs[i, 1])
        x, _, joints, seq = run.path(pores, field, 0)
        V = np.array([run.path(pores, field, k)[1] for k in range(nT)])
        data.append(V)
        ax.plot(x, V[0], color="0.7", lw=1, label="t = 0")
        ln, = ax.plot(x, V[0], lw=2, color="#1f4e9c"); lines.append(ln)
        for j, nd in zip(joints, seq[1:-1]):
            ax.axvline(j, color="k", lw=0.8, ls=":")
            ax.text(j, 1.0, f" node {nd}", transform=ax.get_xaxis_transform(), fontsize=7, va="bottom")
        lo, hi = np.nanmin(V), np.nanmax(V)
        ax.set_ylim(min(0, lo) - 0.04 * (hi - lo), hi + 0.06 * (hi - lo))
        ax.set_xlim(x[0], x[-1])
        ax.set_ylabel(_label(field))
        ax.set_title("path " + " → ".join(f"p{p}" for p in pores) + "  (nodes " + "-".join(map(str, seq)) + ")",
                     fontsize=9, pad=12)
        if i == len(paths) - 1:
            ax.set_xlabel("x along path")
    sup = fig.suptitle("", fontsize=11)

    def upd(k):
        lc.set_array(np.concatenate([run.field(k, p, field) for p in range(run.n_pores)]))
        for ln, V in zip(lines, data):
            ln.set_ydata(V[k])
        sup.set_text(f"{title or ''}   t = {run.times[k]:.4f}".strip())
        return [lc, sup] + lines

    os.makedirs(out, exist_ok=True)
    out_paths = save_animation(fig, upd, nT, os.path.join(out, f"{prefix}_{field}"), gif=gif, mp4=mp4, fps=fps)
    plt.close(fig)
    return out_paths
