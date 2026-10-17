"""Plot functions, one per case. All write PNGs (and optionally a GIF) into `out`."""
import os

import numpy as np

from .io import read_current_flow, read_log, read_network, read_snapshots


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    return plt, animation


# ---------------- single pore: shock propagation ----------------
def front_position(x, c, frac=0.5):
    """First x (from the inlet) where C_bar drops below frac * C_bar(inlet)."""
    below = np.nonzero(c < frac * c[0])[0]
    return x[below[0]] if below.size else np.nan


def single_pore(run_dir, network, dt, write_every, out="docs/figures", prefix="single_pore",
                tmax=np.inf, log=None, pore=0, num_pores=1, gif=True, mp4=False,
                spacetime_cmap="hot"):
    """Space-time map, profiles, front/current plot and GIF for one pore. Returns summary dict."""
    plt, animation = _mpl()
    times, frames = read_snapshots(run_dir, dt, write_every)
    keep = times <= tmax + 1e-12
    times = times[keep]
    C = np.array([f[pore]["C_bar"] for f, k in zip(frames, keep) if k])
    pore = read_network(network, num_pores)["pores"][pore]
    L, ns = pore["length"], C.shape[1]
    x = (np.arange(ns) + 0.5) / ns * L
    xs = np.array([front_position(x, c) for c in C])
    os.makedirs(out, exist_ok=True)
    pre = os.path.join(out, prefix)
    cmap = "viridis"
    tl = "t (nondimensional)"

    # 1) space-time map: x across, t up (linear), black = depleted, white/yellow = bulk,
    #    dashed line = shock front. Same layout as the space-time figure in
    #    Alizadeh, PhD dissertation, Stanford University, 2017 (Fig. 4.2b).
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    m = ax.pcolormesh(x, times, C, shading="auto", cmap=spacetime_cmap, vmin=0,
                      vmax=max(1.0, C.max()), rasterized=True)
    ax.plot(xs, times, "k--", lw=1.2)
    ax.set_xlim(0, L); ax.set_ylim(0, times[-1])
    ax.set_xlabel("x"); ax.set_ylabel("t")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(True)
    cb = fig.colorbar(m, ax=ax, pad=0.03)
    cb.ax.set_title(r"$\bar C$", fontsize=12)
    fig.tight_layout(); fig.savefig(pre + "_spacetime.png", dpi=140); plt.close(fig)

    # 2) waterfall of profiles
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    sel = np.unique(np.r_[0, np.geomspace(1, len(times) - 1, 14).astype(int)])
    colors = plt.get_cmap("plasma")(np.linspace(0, 0.9, len(sel)))
    for c_, k in zip(colors, sel):
        ax.plot(x, C[k], color=c_, lw=1.6, label=f"t = {times[k]:.4g}")
    ax.set_xlabel("x along pore"); ax.set_ylabel(r"$\bar C$")
    ax.set_title(r"$\bar C$ profiles: a sharp front sweeps from the dead end toward the inlet")
    ax.legend(ncol=1, fontsize=8, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout(); fig.savefig(pre + "_profiles.png", dpi=160); plt.close(fig)

    # 3) front position + total current
    fig, ax1 = plt.subplots(figsize=(7.5, 4.2))
    ax1.plot(times, xs, "C0-", lw=2)
    ax1.set_xlabel(tl); ax1.set_ylabel("shock front position x", color="C0")
    ax1.set_ylim(0, L)
    if log:
        tl_, itot = read_log(log)
        k = tl_ <= tmax + 1e-12
        ax2 = ax1.twinx(); ax2.spines["right"].set_visible(True)
        ax2.plot(tl_[k], itot[k], "C3-", lw=1.5)
        ax2.set_ylabel(r"$I_{\rm tot}$", color="C3")
    ax1.set_title("Front position and total current")
    fig.tight_layout(); fig.savefig(pre + "_front_current.png", dpi=160); plt.close(fig)

    # 4) animated GIF
    if gif or mp4:
        fig, ax = plt.subplots(figsize=(6.4, 3.8))
        line, = ax.plot(x, C[0], color="#1f77b4", lw=2.2)
        fill = [ax.fill_between(x, 0, C[0], color="#1f77b4", alpha=0.18)]
        vl = ax.axvline(xs[0], color="0.4", ls="--", lw=1)
        ax.set_xlim(0, L); ax.set_ylim(0, 1.08 * C.max())
        ax.set_xlabel("x along pore (inlet → dead end)"); ax.set_ylabel(r"$\bar C$")
        title = ax.set_title("")

        def upd(k):
            line.set_ydata(C[k])
            fill[0].remove(); fill[0] = ax.fill_between(x, 0, C[k], color="#1f77b4", alpha=0.18)
            vl.set_xdata([xs[k], xs[k]])
            title.set_text(f"Deionization shock, t = {times[k]:.4f}")
            return line, vl, title

        fig.tight_layout()
        from .plots import save_animation
        save_animation(fig, upd, len(times), pre + "_shock", gif=gif, mp4=mp4, fps=16)
        plt.close(fig)

    # summary for the report
    ok = ~np.isnan(xs)
    print(f"frames={len(times)}  t=[{times[0]}, {times[-1]}]  L={L}  cells={ns}")
    for k in sel:
        print(f"  t={times[k]:.4g}  front x={xs[k]:.4f}  C_min={C[k].min():.3g}")
    v = np.gradient(xs[ok], times[ok])
    print(f"front speed dx/dt: initial {v[1]:.1f}, mid {v[len(v)//4]:.1f}, late {v[-1]:.2f}")
    return dict(times=times, front=xs, C=C, x=x)




# ---------------- H-junction: shocks into the branches ----------------
def _channel(frames_k, pores, geom):
    """Concatenate C_bar of a chain of pores; returns x (from the chain start), C, junction x list."""
    xs, cs, joints, x0 = [], [], [], 0.0
    for p in pores:
        L, n = geom[p]["length"], len(frames_k[p]["C_bar"])
        xs.append(x0 + (np.arange(n) + 0.5) / n * L)
        cs.append(frames_k[p]["C_bar"])
        x0 += L
        joints.append(x0)
    return np.concatenate(xs), np.concatenate(cs), joints[:-1]


def _crossing(x, c, start, direction, level):
    """Distance from x=start (walking in direction +1/-1) to where c first crosses `level`
    coming back from the depleted (c < level) side. NaN if c at start is not below level."""
    if direction > 0:
        idx = np.nonzero(x >= start)[0]
    else:
        idx = np.nonzero(x <= start)[0][::-1]
    if idx.size == 0 or c[idx[0]] >= level:
        return np.nan
    above = np.nonzero(c[idx] >= level)[0]
    end = x[idx[above[0]]] if above.size else x[idx[-1]]
    return abs(end - start)


def hjunction(run_dir, network, dt, write_every, out="docs/figures", prefix="hjunction",
              top=(0, 1), bottom=(3, 4), bridge=2, tmax=np.inf, log=None, gif=True,
              snapshot_times=None, mp4=False):
    """H-junction figures: space-time maps of both channels, profiles, network schematic
    coloured by C_bar, depletion/enrichment front distances, currents, GIF."""
    plt, animation = _mpl()
    from matplotlib.collections import LineCollection
    from matplotlib.colors import TwoSlopeNorm

    net = read_network(network)
    geom = net["pores"]
    times, frames = read_snapshots(run_dir, dt, write_every)
    keep = times <= tmax + 1e-12
    times, frames = times[keep], [f for f, k in zip(frames, keep) if k]
    os.makedirs(out, exist_ok=True)
    pre = os.path.join(out, prefix)
    tl = "t (nondimensional)"

    xt, _, jt = _channel(frames[0], top, geom)
    xb, _, jb = _channel(frames[0], bottom, geom)
    Ct = np.array([_channel(f, top, geom)[1] for f in frames])
    Cb = np.array([_channel(f, bottom, geom)[1] for f in frames])
    c_bulk = float(np.median(Ct[0]))
    cmax = max(Ct.max(), Cb.max())
    norm = TwoSlopeNorm(vmin=0.0, vcenter=c_bulk, vmax=cmax)
    cmap = "RdBu_r"
    rlab = lambda pores: " → ".join(f"pore {p}" for p in pores)

    # 1) space-time maps
    fig, axs = plt.subplots(2, 1, figsize=(8.5, 6.6), sharey=True)
    for ax, x, C, j, name, pores in [(axs[0], xt, Ct, jt, "top channel (cation-depleting side)", top),
                                    (axs[1], xb, Cb, jb, "bottom channel (enriching side)", bottom)]:
        m = ax.pcolormesh(x, times, C, shading="auto", cmap=cmap, norm=norm, rasterized=True)
        for jj in j:
            ax.axvline(jj, color="k", lw=0.8, ls=":")
            ax.text(jj, times[-1], " junction", va="top", fontsize=8)
        ax.set_title(f"{name}: {rlab(pores)}", fontsize=10)
        ax.set_ylabel(tl)
    axs[1].set_xlabel("x along channel")
    fig.colorbar(m, ax=axs, label=r"$\bar C$ (white = bulk)")
    fig.suptitle(r"H-junction: $\bar C(x,t)$ in the two micro-channels")
    fig.savefig(pre + "_spacetime.png", dpi=160, bbox_inches="tight"); plt.close(fig)

    # 2) profiles at selected times
    sel = np.unique(np.r_[0, np.geomspace(1, len(times) - 1, 8).astype(int)])
    colors = plt.get_cmap("plasma")(np.linspace(0, 0.9, len(sel)))
    fig, axs = plt.subplots(2, 1, figsize=(8.5, 6.2))
    for ax, x, C, j, name in [(axs[0], xt, Ct, jt, "top channel"), (axs[1], xb, Cb, jb, "bottom channel")]:
        for c_, k in zip(colors, sel):
            ax.plot(x, C[k], color=c_, lw=1.4, label=f"t = {times[k]:.3g}")
        for jj in j:
            ax.axvline(jj, color="k", lw=0.8, ls=":")
        ax.set_ylabel(r"$\bar C$"); ax.set_title(name, fontsize=10)
    axs[1].set_xlabel("x along channel (dotted: junction with the nano-bridge)")
    axs[0].legend(ncol=1, fontsize=8, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.0))
    fig.tight_layout(); fig.savefig(pre + "_profiles.png", dpi=160, bbox_inches="tight"); plt.close(fig)

    # 3) schematic coloured by C_bar at selected times
    if snapshot_times is None:
        snapshot_times = [times[0], times[len(times) // 10], times[len(times) // 3], times[-1]]
    kk = [int(np.argmin(abs(times - t))) for t in snapshot_times]
    Lt0, Lb0 = geom[top[0]]["length"], geom[bottom[0]]["length"]
    yb = -0.6

    def segs(f):
        S, V = [], []
        for pores, y, xoff in [(top, 0.0, 0.0), (bottom, yb, Lt0 - Lb0)]:
            x, c, _ = _channel(f, pores, geom)
            x = x + xoff
            edges = np.r_[x[0] - (x[1] - x[0]) / 2, (x[1:] + x[:-1]) / 2, x[-1] + (x[-1] - x[-2]) / 2]
            S += [[(edges[i], y), (edges[i + 1], y)] for i in range(len(x))]
            V += list(c)
        c = f[bridge]["C_bar"]
        ys = np.linspace(0, yb, len(c) + 1)
        S += [[(Lt0, ys[i]), (Lt0, ys[i + 1])] for i in range(len(c))]
        V += list(c)
        return S, np.array(V)

    fig, axs = plt.subplots(2, 2, figsize=(10, 5.2))
    for ax, k in zip(axs.flat, kk):
        S, V = segs(frames[k])
        lc = LineCollection(S, cmap=cmap, norm=norm, linewidths=7, capstyle="butt")
        lc.set_array(V); ax.add_collection(lc)
        res = net["reservoirs"]
        ends_xy = {top[0]: None}
        for pores, y, xoff in [(top, 0.0, 0.0), (bottom, yb, Lt0 - Lb0)]:
            xx = xoff
            nodes = [geom[pores[0]]["ends"][0]] + [geom[p]["ends"][1] for p in pores]
            pos = [xx]
            for p in pores:
                xx += geom[p]["length"]; pos.append(xx)
            for nd, px in zip(nodes, pos):
                fixed = res and res[nd]["potential_known"]
                ax.plot(px, y, "s" if fixed else "o", ms=8 if fixed else 6,
                        color="k" if fixed else "w", mec="k", zorder=3)
                if fixed and k == kk[0]:
                    ax.annotate(f"V={res[nd]['potential']:g}", (px, y), textcoords="offset points",
                                xytext=(0, 9), ha="center", fontsize=8)
        ax.set_xlim(Lt0 - Lb0 - 0.1, Lt0 + Lb0 + 0.1); ax.set_ylim(yb - 0.25, 0.3)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"t = {times[k]:.3g}", fontsize=10)
    fig.colorbar(lc, ax=axs, shrink=0.8, label=r"$\bar C$")
    fig.suptitle("H-junction network coloured by $\\bar C$ (nano-bridge drawn long, not to scale;\n"
                 "squares = fixed-potential reservoirs, circles = internal junctions)", fontsize=10)
    fig.savefig(pre + "_schematic.png", dpi=160, bbox_inches="tight"); plt.close(fig)

    # 4) front distances from the junctions
    level = 0.5 * c_bulk
    d_up = np.array([_crossing(xt, c, jt[0], -1, level) for c in Ct])
    d_dn = np.array([_crossing(xt, c, jt[0], +1, level) for c in Ct])
    enr = 0.5 * (c_bulk + Cb.max(axis=1))
    e_up = np.array([abs(jb[0] - xb[np.nonzero((xb <= jb[0]) & (c < max(e, 1.05 * c_bulk)))[0][-1]])
                     if (c > 1.05 * c_bulk).any() else np.nan for c, e in zip(Cb, enr)])
    fig, ax = plt.subplots(figsize=(8, 4.4))
    ax.plot(times, d_up, lw=2, label=f"depletion front, upstream into pore {top[0]}")
    ax.plot(times, d_dn, lw=2, label=f"depletion front, downstream into pore {top[1]}")
    ax.plot(times, e_up, lw=2, ls="--", label=f"enriched zone half-width, pore {bottom[0]}")
    ax.set_xlabel(tl); ax.set_ylabel("distance from junction")
    ax.set_title(r"Fronts leaving the junctions ($\bar C = \frac{1}{2}\bar C_{\rm bulk}$ for depletion)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(pre + "_fronts.png", dpi=160); plt.close(fig)

    # 5) currents
    cf_path = os.path.join(run_dir, "current_flow_area.bin")
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    if log:
        t_, itot = read_log(log)
        k_ = t_ <= tmax + 1e-12
        axs[0].plot(t_[k_], itot[k_], "C3-", lw=1.5)
        axs[0].set_xlabel(tl); axs[0].set_ylabel(r"$I_{\rm tot}$"); axs[0].set_title("total current", fontsize=10)
    if os.path.exists(cf_path):
        cf = read_current_flow(cf_path)
        idx = np.arange(len(cf["I"]))
        axs[1].bar(idx - 0.2, cf["I"], 0.4, label="current I")
        axs[1].bar(idx + 0.2, cf["Q"], 0.4, label="flow rate Q")
        axs[1].axhline(0, color="k", lw=0.6)
        axs[1].set_xticks(idx, [f"pore {i}" for i in idx])
        axs[1].set_title("per-pore current and flow at the final time", fontsize=10)
        axs[1].legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(pre + "_currents.png", dpi=160); plt.close(fig)

    # 6) GIF
    if gif or mp4:
        fig, axs = plt.subplots(2, 1, figsize=(6.8, 4.8), sharex=False)
        lines = []
        for ax, x, C, j, name in [(axs[0], xt, Ct, jt, "top channel"), (axs[1], xb, Cb, jb, "bottom channel")]:
            ln, = ax.plot(x, C[0], lw=2, color="#1f77b4"); lines.append(ln)
            for jj in j:
                ax.axvline(jj, color="k", lw=0.8, ls=":")
            ax.set_ylim(min(0, C.min()) - 0.05, 1.08 * C.max()); ax.set_ylabel(r"$\bar C$")
            ax.set_title(name, fontsize=9)
        title = fig.suptitle("")
        fig.tight_layout()

        def upd(k):
            lines[0].set_ydata(Ct[k]); lines[1].set_ydata(Cb[k])
            title.set_text(f"H-junction, t = {times[k]:.4f}")
            return lines + [title]

        from .plots import save_animation
        save_animation(fig, upd, len(times), pre + "_shock", gif=gif, mp4=mp4, fps=12)
        plt.close(fig)

    print(f"frames={len(times)} t=[{times[0]}, {times[-1]}]  C_bulk={c_bulk:.3f}  C_max(bottom)={Cb.max():.3f}")
    for k in sel:
        print(f"  t={times[k]:.4g}  upstream front={d_up[k]:.3f}  downstream front={d_dn[k]:.3f}  "
              f"enriched half-width={e_up[k]:.3f}  junction C bottom={Cb[k][np.argmin(abs(xb - jb[0]))]:.2f}")
    if os.path.exists(cf_path):
        cf = read_current_flow(cf_path)
        print("  final per-pore I:", np.round(cf["I"], 3), " Q:", np.round(cf["Q"], 4))
    return dict(times=times, d_up=d_up, d_dn=d_dn)
