"""VTK export (ParaView, Blender importers, PyVista) and an off-screen 3D tube renderer.

export_vtk writes one PolyData file per snapshot (<prefix>_NNNN.vtp, XML, ASCII) plus <prefix>.pvd,
a time series ParaView opens directly. Geometry: every pore is a straight polyline between its two
nodes (coordinates from x/y in network.toml, Network.in geometry rows, a CSV, or the spring layout);
each computational cell is one line segment. Cell data per segment: C_bar, C0, P, Mu, pore_id.
No external dependency (plain XML), so it scales to large random networks.

render3d (needs pyvista) draws the pores as tubes coloured by a field and writes PNG frames and/or
an MP4, off-screen.
"""
import os
from xml.sax.saxutils import quoteattr

import numpy as np

from .layout import node_positions

FIELDS = ("C_bar", "C0", "P", "Mu")


def geometry(run, layout=None):
    """Points (N,3), segments (M,2) and pore id per segment; segment order = pores then cells."""
    pos, src = node_positions(run, layout)
    pts, segs, pid = [], [], []
    for p, pore in enumerate(run.net["pores"]):
        a, b = np.r_[pos[pore["ends"][0]], 0.0], np.r_[pos[pore["ends"][1]], 0.0]
        n = len(run.field(0, p, "C_bar"))
        s = np.linspace(0.0, 1.0, n + 1)[:, None]
        base = len(pts)
        pts.extend(a + s * (b - a))
        segs.extend([(base + i, base + i + 1) for i in range(n)])
        pid.extend([p] * n)
    return np.asarray(pts), np.asarray(segs, dtype=np.int64), np.asarray(pid, dtype=np.int32), src


def _cell_fields(run, k, fields):
    out = {}
    for f in fields:
        try:
            out[f] = np.concatenate([np.asarray(run.field(k, p, f), float) for p in range(run.n_pores)])
        except Exception:  # noqa: BLE001  field not in this snapshot
            pass
    return out


def _arr(name, a, ncomp=1, typ="Float64"):
    fmt = "%d" if typ.startswith("Int") else "%.9g"
    body = " ".join(fmt % v for v in np.ravel(a))
    return (f'<DataArray type="{typ}" Name={quoteattr(name)} NumberOfComponents="{ncomp}" format="ascii">'
            f"{body}</DataArray>\n")


def write_vtp(path, pts, segs, cell_data):
    with open(path, "w") as fh:
        fh.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n'
                 "<PolyData>\n")
        fh.write(f'<Piece NumberOfPoints="{len(pts)}" NumberOfVerts="0" NumberOfLines="{len(segs)}" '
                 'NumberOfStrips="0" NumberOfPolys="0">\n')
        fh.write("<Points>\n" + _arr("Points", pts, 3) + "</Points>\n")
        fh.write("<CellData Scalars=\"C_bar\">\n")
        for name, a in cell_data.items():
            fh.write(_arr(name, a, 1, "Int32" if a.dtype.kind == "i" else "Float64"))
        fh.write("</CellData>\n<Lines>\n")
        fh.write(_arr("connectivity", segs, 1, "Int64"))
        fh.write(_arr("offsets", 2 * np.arange(1, len(segs) + 1), 1, "Int64"))
        fh.write("</Lines>\n</Piece>\n</PolyData>\n</VTKFile>\n")


def export_vtk(run, out, prefix="network", layout=None, fields=FIELDS, every=1):
    """Write <out>/<prefix>_NNNN.vtp for every `every`-th snapshot and <out>/<prefix>.pvd."""
    os.makedirs(out, exist_ok=True)
    pts, segs, pid, src = geometry(run, layout)
    entries = []
    for i, k in enumerate(range(0, len(run.times), every)):
        cd = _cell_fields(run, k, fields)
        cd["pore_id"] = pid
        name = f"{prefix}_{i:04d}.vtp"
        write_vtp(os.path.join(out, name), pts, segs, cd)
        entries.append((float(run.times[k]), name))
    pvd = os.path.join(out, f"{prefix}.pvd")
    with open(pvd, "w") as fh:
        fh.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n<Collection>\n')
        for t, name in entries:
            fh.write(f'<DataSet timestep="{t:.9g}" part="0" file="{name}"/>\n')
        fh.write("</Collection>\n</VTKFile>\n")
    return pvd, len(entries), src


def render3d(run, out, prefix="render3d", field="C_bar", layout=None, mp4=True, png_times=None,
             radius=None, fps=12, cmap="RdBu_r", clim=None, window=(1280, 720), title=None,
             log_scale=None):
    """Pores as tubes coloured by `field` (cell data interpolated to points), rendered off-screen.
    Writes <prefix>_<field>.mp4 (every snapshot) and PNG frames at `png_times`. Needs pyvista."""
    import pyvista as pv
    pv.OFF_SCREEN = True
    pts, segs, pid, src = geometry(run, layout)
    lines = np.hstack([np.full((len(segs), 1), 2), segs]).ravel()
    mesh = pv.PolyData(pts, lines=lines)
    span = float(np.ptp(pts[:, :2], axis=0).max()) or 1.0
    radius = radius or 0.012 * span
    allv = np.concatenate([_cell_fields(run, k, [field])[field] for k in range(len(run.times))])
    if clim is None:
        lo, hi = np.percentile(allv, [0.5, 99.5])
        clim = (float(lo), float(hi))
    # depletion (~0) together with strong enrichment (>> 1): a log colour scale shows both
    log_scale = bool(log_scale) if log_scale is not None else bool(clim[0] > 0 and clim[1] > 3.0 and clim[1] / clim[0] > 30)
    if log_scale:
        clim = (max(clim[0], 1e-3), clim[1])
    os.makedirs(out, exist_ok=True)
    base = os.path.join(out, f"{prefix}_{field}")
    written = []

    def frame_mesh(k):
        m = mesh.copy()
        m.cell_data[field] = _cell_fields(run, k, [field])[field]
        return m.cell_data_to_point_data().tube(radius=radius, n_sides=16)

    def setup(pl, k):
        tube = frame_mesh(k)
        pl.add_mesh(tube, scalars=field, cmap=cmap, clim=clim, smooth_shading=True, log_scale=log_scale,
                    scalar_bar_args=dict(title=field, vertical=True, position_x=0.88, position_y=0.2))
        res = run.net.get("reservoirs") or []
        pos, _ = node_positions(run, layout)
        for nd in range(run.n_nodes):
            fixed = bool(res) and nd < len(res) and res[nd]["potential_known"]
            c = (float(pos[nd][0]), float(pos[nd][1]), 0.0)
            glyph = pv.Box(bounds=(c[0] - 1.8 * radius, c[0] + 1.8 * radius, c[1] - 1.8 * radius,
                                   c[1] + 1.8 * radius, -1.8 * radius, 1.8 * radius)) \
                if fixed else pv.Sphere(radius=1.6 * radius, center=c)
            pl.add_mesh(glyph, color="black" if fixed else "white", smooth_shading=not fixed)
        pl.view_xy(); pl.camera.elevation = -35; pl.camera.zoom(1.1); pl.reset_camera_clipping_range()
        pl.set_background("white")
        return tube

    if mp4:
        pl = pv.Plotter(off_screen=True, window_size=window)
        tube = setup(pl, 0)
        txt = pl.add_text("", font_size=11, color="black")
        pl.open_movie(base + ".mp4", framerate=fps, quality=7)
        for k in range(len(run.times)):
            tube.copy_from(frame_mesh(k))   # same dataset object -> the actor updates in place
            txt.SetText(2, f"{title or ''}  t = {run.times[k]:.4f}".strip())
            pl.write_frame()
        pl.close()
        written.append(base + ".mp4")
    for t in (png_times or [run.times[-1]]):
        k = int(np.argmin(abs(run.times - t)))
        pl = pv.Plotter(off_screen=True, window_size=window)
        setup(pl, k)
        pl.add_text(f"{title or ''}  t = {run.times[k]:.4g}".strip(), font_size=11, color="black")
        p = f"{base}_t{run.times[k]:.4g}.png"
        pl.screenshot(p); pl.close()
        written.append(p)
    return written
