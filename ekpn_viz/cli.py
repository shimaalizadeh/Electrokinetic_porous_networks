"""ekpn-viz: plots for any ek_solver network.

Generic commands (work for any network):
  ekpn-viz network RUN_DIR ...  schematic coloured by a field at chosen times
  ekpn-viz path    RUN_DIR ... --pores 0,2,1   field along a chain of pores (+ fronts, GIF/MP4)
  ekpn-viz export-vtk RUN_DIR ...  .vtp per snapshot + .pvd (ParaView / Blender import)
  ekpn-viz render3d   RUN_DIR ...  tubes coloured by C_bar, off-screen MP4 + PNG (pyvista)
  ekpn-viz animate RUN_DIR ... --paths "0,1;3,4"  network schematic + paths, MP4 and small GIF
  ekpn-viz pores   RUN_DIR ...  space-time map of every pore
  ekpn-viz flow    RUN_DIR ...  U map with flow arrows, sign reversals, circulating loops
  ekpn-viz current RUN_DIR ...  I_tot(t) and final per-pore current/flow, Kirchhoff check
Presets (fixed figure sets used in docs/figures): single-pore, hjunction.

All commands need --network (Network.in or network.toml used for the run), --dt and
--write-every (to turn snapshot numbers into times).
"""
import argparse

import numpy as np

from . import plots, presets
from .io import Run


def _ints(s):
    return [int(v) for v in s.replace(" ", "").split(",") if v != ""]


def _floats(s):
    return [float(v) for v in s.replace(" ", "").split(",") if v != ""]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ekpn-viz", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(s, prefix, field=True):
        s.add_argument("run_dir", help="folder with Tecplot_*.bin (the solver's working directory)")
        s.add_argument("--network", required=True, help="Network.in or network.toml used for the run")
        s.add_argument("--dt", type=float, required=True)
        s.add_argument("--write-every", type=int, required=True)
        s.add_argument("--num-pores", type=int, help="legacy Network.in only (default: from snapshots)")
        s.add_argument("--num-reservoirs", type=int, help="legacy Network.in only, for labels/layout")
        s.add_argument("--tmax", type=float, default=np.inf, help="only use frames with t <= tmax")
        s.add_argument("--log", help="solver stdout log, for I_tot(t)")
        s.add_argument("--out", default="docs/figures")
        s.add_argument("--prefix", default=prefix)
        if field:
            s.add_argument("--field", default="C_bar", choices=["C_bar", "C0", "P", "Mu", "U"])

    s = sub.add_parser("network", help="schematic coloured by a field")
    common(s, "network")
    s.add_argument("--times", type=_floats, help="comma-separated times (default: 4 spread over the run)")
    s.add_argument("--shape", choices=["line", "slit"], default="line",
                   help="pores as coloured lines or as flat slit channels (oblique slab look)")
    s.add_argument("--title")
    s.add_argument("--layout", help="CSV node,x,y or 'spring' (default: x/y in network.toml, "
                                    "Network.in geometry rows, else spring layout)")
    s = sub.add_parser("path", help="field along a chain of pores")
    common(s, "path")
    s.add_argument("--pores", type=_ints, required=True, help="e.g. 0,2,1 (must share nodes in order)")
    s.add_argument("--threshold", type=float, help="front level (default: half the initial C_bar)")
    s.add_argument("--no-gif", action="store_true")
    s.add_argument("--mp4", action="store_true", help="also write an MP4 (needs ffmpeg)")
    s.add_argument("--log-time", action="store_true", help="log time axis in the space-time map and front plot")
    s = sub.add_parser("pores", help="space-time map of every pore")
    common(s, "pores")
    s.add_argument("--pores", type=_ints, help="subset (default: all)")
    s = sub.add_parser("flow", help="U on the network with direction arrows, reversals, circulating loops")
    common(s, "flow", field=False)
    s.add_argument("--times", type=_floats)
    s.add_argument("--layout")
    s = sub.add_parser("current", help="I_tot(t) and per-pore current/flow")
    common(s, "currents", field=False)

    s = sub.add_parser("animate", help="network schematic + field along paths, MP4 and small GIF")
    common(s, "anim")
    s.add_argument("--paths", required=True, help='pore chains separated by ";", e.g. "0,1;3,4"')
    s.add_argument("--shape", choices=["line", "slit"], default="line")
    s.add_argument("--layout")
    s.add_argument("--fps", type=int, default=12)
    s.add_argument("--title")
    s.add_argument("--no-gif", action="store_true")
    s.add_argument("--no-mp4", action="store_true")
    s = sub.add_parser("export-vtk", help="VTK PolyData (.vtp) per snapshot + .pvd time series (ParaView/Blender)")
    common(s, "network", field=False)
    s.add_argument("--layout"); s.add_argument("--every", type=int, default=1, help="export every N-th snapshot")
    s = sub.add_parser("render3d", help="pores as tubes coloured by a field, off-screen MP4/PNG (needs pyvista)")
    common(s, "render3d")
    s.add_argument("--layout"); s.add_argument("--times", type=_floats, help="PNG frames at these times")
    s.add_argument("--no-mp4", action="store_true"); s.add_argument("--fps", type=int, default=12)
    s.add_argument("--title")
    s = sub.add_parser("single-pore", help="preset: figures in docs/figures/single_pore_*")
    common(s, "single_pore", field=False)
    s.add_argument("--pore", type=int, default=0)
    s.add_argument("--no-gif", action="store_true")
    s.add_argument("--mp4", action="store_true", help="also write an MP4 (needs ffmpeg)")
    s = sub.add_parser("hjunction", help="preset: figures in docs/figures/hjunction_*")
    common(s, "hjunction", field=False)
    s.add_argument("--top", type=_ints, default=[0, 1])
    s.add_argument("--bottom", type=_ints, default=[3, 4])
    s.add_argument("--bridge", type=int, default=2)
    s.add_argument("--no-gif", action="store_true")
    s.add_argument("--mp4", action="store_true", help="also write an MP4 (needs ffmpeg)")

    a = ap.parse_args(argv)
    if a.cmd == "single-pore":
        presets.single_pore(a.run_dir, a.network, a.dt, a.write_every, out=a.out, prefix=a.prefix,
                            tmax=a.tmax, log=a.log, pore=a.pore, num_pores=a.num_pores or 1,
                            gif=not a.no_gif, mp4=a.mp4)
        return
    if a.cmd == "hjunction":
        presets.hjunction(a.run_dir, a.network, a.dt, a.write_every, out=a.out, prefix=a.prefix,
                          top=a.top, bottom=a.bottom, bridge=a.bridge, tmax=a.tmax, log=a.log,
                          gif=not a.no_gif, mp4=a.mp4)
        return
    run = Run(a.run_dir, a.network, a.dt, a.write_every, a.tmax, a.num_pores, a.num_reservoirs, a.log)
    if a.cmd == "network":
        print(plots.network_map(run, a.field, a.times, a.out, a.prefix, a.layout, shape=a.shape, title=a.title))
    elif a.cmd == "path":
        r = plots.path_plots(run, a.pores, a.field, a.threshold, a.out, a.prefix, gif=not a.no_gif, mp4=a.mp4,
                             log_time=a.log_time)
        print(f"path nodes {r['nodes']}, length {r['x'][-1]:.4g}, {len(r['fronts'])} front points")
    elif a.cmd == "pores":
        print(plots.pores_grid(run, a.field, a.out, a.prefix, a.pores))
    elif a.cmd == "flow":
        print(plots.network_map(run, "U", a.times, a.out, a.prefix, a.layout))
        path, rep = plots.flow_report(run, a.out, a.prefix)
        print(path)
        for line in rep:
            print("  " + line)
    elif a.cmd == "export-vtk":
        from . import vtk
        pvd, n, src = vtk.export_vtk(run, a.out, a.prefix, a.layout, every=a.every)
        print(f"{pvd}: {n} snapshots (layout: {src})")
    elif a.cmd == "render3d":
        from . import vtk
        for pth in vtk.render3d(run, a.out, a.prefix, a.field, a.layout, mp4=not a.no_mp4,
                                png_times=a.times, fps=a.fps, title=a.title):
            print(pth)
    elif a.cmd == "animate":
        paths = [_ints(c) for c in a.paths.split(";")]
        for pth in plots.animate(run, paths, a.field, a.out, a.prefix, a.layout, gif=not a.no_gif,
                                 mp4=not a.no_mp4, fps=a.fps, title=a.title, shape=a.shape):
            print(pth)
    elif a.cmd == "current":
        path, cf = plots.currents(run, a.out, a.prefix)
        print(path)
        if cf is not None:
            print("per-pore I:", np.round(cf["I"], 4), " Q:", np.round(cf["Q"], 4))
            print("node balance (in - out):", np.round(plots.node_balance(run, cf), 4))


if __name__ == "__main__":
    main()
