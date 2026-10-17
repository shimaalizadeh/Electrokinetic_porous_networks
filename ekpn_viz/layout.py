"""Node layout for network schematics: given coordinates, or a deterministic spring layout."""
import numpy as np


def spring_layout(n_nodes, edges, lengths=None, iters=500, seed=0):
    """Fruchterman-Reingold style layout in plain numpy. edges: list of (a, b)."""
    rng = np.random.default_rng(seed)
    ang = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    pos = np.c_[np.cos(ang), np.sin(ang)] + 0.01 * rng.standard_normal((n_nodes, 2))
    if n_nodes <= 2:
        return np.array([[0.0, 0.0], [1.0, 0.0]])[:n_nodes]
    k = 1.0 / np.sqrt(n_nodes)
    t = 0.1
    E = np.array(edges)
    for _ in range(iters):
        d = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(d, axis=-1) + 1e-9
        disp = (d / dist[..., None] * (k * k / dist)[..., None]).sum(1)       # repulsion
        dv = pos[E[:, 0]] - pos[E[:, 1]]
        dl = np.linalg.norm(dv, axis=1)[:, None] + 1e-9
        f = dv / dl * dl ** 2 / k                                                # attraction
        np.add.at(disp, E[:, 0], -f)
        np.add.at(disp, E[:, 1], f)
        ln = np.linalg.norm(disp, axis=1)[:, None] + 1e-9
        pos += disp / ln * np.minimum(ln, t)
        t *= 0.99
    pos -= pos.mean(0)
    return pos / np.abs(pos).max()


def read_layout_csv(path):
    """CSV 'node,x,y' (header optional) -> {node: (x, y)}."""
    out = {}
    for line in open(path):
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 3 or not parts[0].lstrip("-").isdigit():
            continue
        out[int(parts[0])] = (float(parts[1]), float(parts[2]))
    return out


def node_positions(run, layout=None):
    """layout: None (auto), 'spring', a CSV path, or a dict. Returns array (n_nodes, 2), source str."""
    n = run.n_nodes
    edges = [tuple(p["ends"]) for p in run.net["pores"]]
    coords = None
    src = "spring layout"
    if isinstance(layout, dict):
        coords, src = layout, "given"
    elif layout and layout != "spring":
        coords, src = read_layout_csv(layout), layout
    elif layout is None and run.net.get("coords"):
        coords, src = run.net["coords"], "network.toml x/y"
    elif layout is None and not str(run.network_path).endswith(".toml"):
        from .io import legacy_geometry
        coords = legacy_geometry(run.network_path, run.n_pores, n)
        src = "Network.in geometry rows"
    if coords and all(i in coords for i in range(n)) and len({coords[i] for i in range(n)}) == n:
        return np.array([coords[i] for i in range(n)], float), src
    return spring_layout(n, edges), "spring layout"
