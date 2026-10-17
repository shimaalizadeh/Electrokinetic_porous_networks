#!/usr/bin/env python3
"""Shock distance d_s(t) in Channel 1 of the H-junction, defined as in Shima's 2018 scripts
(Drive hjunction/shock_propagation.py, shocks_distance.py): take C_bar along Channel 1 + Channel 2,
level = (min + max)/2 over Channel 1 + first half of Channel 2, scan left from the middle of
Channel 2 while C_bar < level; d_s = (x_junction - x_found) * Lref.

    python3 tools/hjunction_ds.py RUN_DIR --dt 1e-5 --write-every 1000 [--lref-um 6594] [--csv out.csv]
"""
import argparse, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ekpn_viz.io import read_snapshots  # noqa: E402


def ds_of(frame, lref_um, L1=1.0):
    c1 = np.asarray(frame[0]["C_bar"]); c2 = np.asarray(frame[1]["C_bar"])
    n1, n2 = len(c1), len(c2)
    c = np.concatenate([c1, c2])
    x = np.concatenate([(np.arange(n1) + 0.5) * L1 / n1, L1 + (np.arange(n2) + 0.5) * L1 / n2])
    idx = n1 + n2 // 2
    level = 0.5 * (c[:idx].min() + c[:idx].max())
    while idx > 0 and c[idx] < level:
        idx -= 1
    return (L1 - x[idx]) * lref_um, float(c1.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir"); ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--write-every", type=int, required=True)
    ap.add_argument("--lref-um", type=float, default=6594.0); ap.add_argument("--csv"); ap.add_argument("--png"); ap.add_argument("--ref-um", type=float)
    a = ap.parse_args()
    times, frames = read_snapshots(a.run_dir, a.dt, a.write_every)
    rows = [(t, *ds_of(f, a.lref_um)) for t, f in zip(times, frames)]
    for t, d, m in rows:
        print(f"{t:10.5f}  ds = {d:9.2f} um   min C_bar(ch1) = {m:.4g}")
    if a.png:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        r = np.array([row for row in rows if row[0] > 0 and row[1] > 0])
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        ax.plot(r[:, 0], r[:, 1], "k-", lw=2, label="shock 1 (Channel 1), this solver")
        if a.ref_um:
            ax.axhline(a.ref_um, color="C3", ls="--", label=f"reference value {a.ref_um:g} um")
        ax.set_xscale("log"); ax.set_xlabel("t (nondimensional)"); ax.set_ylabel("$d_s$ from junction (um)")
        ax.set_title("Deionization shock distance from the cation-selective element", fontsize=10)
        ax.legend(frameon=False); fig.tight_layout(); fig.savefig(a.png, dpi=160)
    if a.csv:
        np.savetxt(a.csv, rows, delimiter=",", header="t,ds_um,minC_ch1", comments="")


if __name__ == "__main__":
    main()
