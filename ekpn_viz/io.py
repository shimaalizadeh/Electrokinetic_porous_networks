"""Readers for ek_solver inputs and outputs.

Output files written by ek_solver into the working directory:
  Tecplot_XXXX.bin  one snapshot per write_every steps (frame k at t = k*write_every*dt;
                    Tecplot_0000 is t=0, the last file is the final state, a repeat of
                    the previous frame when max_time is a multiple of the write interval).
                    Layout (network_write.hpp): int32 n_pores; n_pores x (int32 Ns, int32 nvar=7);
                    then per pore 7 blocks of Ns float64: C_bar, C0, P, Phi(mu), u, v, w.
                    Interior cells only (the two reservoir cells are not written).
  current_flow_area.bin  int32 n; float64 Q[n], I[n], S[n] (final per-pore flow, current, area)
  stdout log        "time= t" lines followed by "I_tot = value"
"""
import glob
import os
import re
import sys

import numpy as np

VARS = ["C_bar", "C0", "P", "Phi", "u", "v", "w"]


# ---------------- readers ----------------
def read_tecplot(path):
    """Return a list (one per pore) of dicts VAR -> array of length Ns."""
    raw = open(path, "rb").read()
    off = 0
    n = int(np.frombuffer(raw, np.int32, 1, off)[0]); off += 4
    hdr = np.frombuffer(raw, np.int32, 2 * n, off).reshape(n, 2); off += 8 * n
    out = []
    for ns, nv in hdr:
        a = np.frombuffer(raw, np.float64, ns * nv, off).reshape(nv, ns); off += 8 * ns * nv
        out.append(dict(zip(VARS, a)))
    if off != len(raw):
        raise ValueError(f"{path}: {len(raw) - off} unexpected trailing bytes")
    return out


def read_snapshots(run_dir, dt, write_every, drop_final_repeat=True):
    """All Tecplot frames, sorted. Returns (times, frames)."""
    files = sorted(glob.glob(os.path.join(run_dir, "Tecplot_*.bin")),
                   key=lambda f: int(re.findall(r"(\d+)\.bin$", f)[0]))
    if not files:
        sys.exit(f"no Tecplot_*.bin in {run_dir}")
    frames = [read_tecplot(f) for f in files]
    idx = np.array([int(re.findall(r"(\d+)\.bin$", f)[0]) for f in files])
    times = idx * write_every * dt
    if drop_final_repeat and len(frames) > 1 and all(
            np.array_equal(a["C_bar"], b["C_bar"]) for a, b in zip(frames[-1], frames[-2])):
        frames, times = frames[:-1], times[:-1]
    return times, frames


def read_current_flow(path):
    raw = open(path, "rb").read()
    n = int(np.frombuffer(raw, np.int32, 1, 0)[0])
    a = np.frombuffer(raw, np.float64, 3 * n, 4).reshape(3, n)
    return dict(Q=a[0], I=a[1], S=a[2])


def read_log(path):
    t, i_tot, cur = [], [], None
    for line in open(path):
        if line.startswith("time="):
            cur = float(line.split()[1])
        elif line.startswith("I_tot") and cur is not None:
            t.append(cur); i_tot.append(float(line.split()[-1])); cur = None
    return np.array(t), np.array(i_tot)


def read_network(path, num_pores=None, num_reservoirs=None):
    """Network geometry from network.toml or legacy Network.in.

    Returns dict(pores=[{length, n_cells, ends}], reservoirs=[{potential, pressure_known,
    potential_known, volume}]). For a legacy file the reservoir list needs num_reservoirs.
    """
    if str(path).endswith(".toml"):
        try:
            import tomllib
        except ModuleNotFoundError:  # Python < 3.11
            import tomli as tomllib
        d = tomllib.load(open(path, "rb"))
        dflt = d.get("defaults", {})
        pores = sorted(d["pores"], key=lambda p: p["id"])
        res = sorted(d.get("reservoirs", []), key=lambda r: r["id"])
        return dict(
            pores=[dict(length=p.get("length", dflt.get("length")),
                        n_cells=p.get("n_cells", dflt.get("n_cells")), ends=list(p["ends"]),
                        area=float(p.get("area", dflt.get("area", 1.0))),
                        type=p.get("type", dflt.get("type", "slit"))) for p in pores],
            reservoirs=[dict(potential=r["potential"], pressure_known=r["pressure_known"],
                             potential_known=r["potential_known"], volume=r.get("volume", 0.0)) for r in res],
            coords={r["id"]: (r["x"], r["y"]) for r in res if "x" in r and "y" in r} or None)
    toks = open(path).read().split()
    n = num_pores or 1
    pores = [dict(length=float(toks[10 * i + 4]), n_cells=int(toks[10 * i + 5]),
                  ends=[int(toks[10 * i + 2]), int(toks[10 * i + 3])]) for i in range(n)]
    res = []
    for j in range(num_reservoirs or 0):
        t = toks[10 * n + 9 * j: 10 * n + 9 * j + 9]
        res.append(dict(potential=float(t[4]), pressure_known=int(t[1]) == 0,
                        potential_known=int(t[2]) == 0, volume=float(t[8])))
    return dict(pores=pores, reservoirs=res)


# ---------------- generic run access ----------------
FIELD_ALIASES = {"C_bar": "C_bar", "C": "C_bar", "C0": "C0", "P": "P", "Mu": "Phi", "Phi": "Phi", "U": "U"}
FIELD_LABELS = {"C_bar": r"$\bar C$", "C0": r"$C_0$", "P": "P", "Mu": r"$\mu$", "Phi": r"$\mu$", "U": "U"}


def legacy_geometry(path, num_pores, num_reservoirs):
    """Node coordinates from the trailing geometry rows some Network.in files carry
    (10 values per row: Ns Nt Nz ds dn dz x0 y0 z0 theta; rows for pores then reservoirs).
    Works for the current and the v0 format. Returns {node: (x, y)} or None."""
    toks = open(path).read().split()
    n_geo = 10 * (num_pores + num_reservoirs)
    base_now = 10 * num_pores + 9 * num_reservoirs + 3 * num_pores
    base_v0 = 1 + 4 * num_pores + 8 * num_reservoirs + 3 * num_pores
    if len(toks) - n_geo not in (base_now, base_v0):
        return None
    geo = toks[len(toks) - n_geo:]
    try:
        rows = [[float(v) for v in geo[10 * i:10 * i + 10]] for i in range(num_pores + num_reservoirs)]
    except ValueError:
        return None
    return {j: (rows[num_pores + j][6], rows[num_pores + j][7]) for j in range(num_reservoirs)}


class Run:
    """One solver run: network + snapshots. Fields per pore: C_bar, C0, P, Mu, U."""

    def __init__(self, run_dir, network, dt, write_every, tmax=np.inf, num_pores=None,
                 num_reservoirs=None, log=None):
        self.run_dir, self.network_path, self.log = run_dir, network, log
        self.times, frames = read_snapshots(run_dir, dt, write_every)
        keep = self.times <= tmax + 1e-12
        self.times = self.times[keep]
        self.frames = [f for f, k in zip(frames, keep) if k]
        n_p = num_pores or len(self.frames[0])
        self.net = read_network(network, n_p, num_reservoirs)
        if len(self.net["pores"]) != len(self.frames[0]):
            raise ValueError(f"network has {len(self.net['pores'])} pores, snapshots have {len(self.frames[0])}")
        self.n_pores = n_p
        self.n_nodes = 1 + max(max(p["ends"]) for p in self.net["pores"])

    def field(self, k, pore, name):
        f = self.frames[k][pore]
        key = FIELD_ALIASES.get(name, name)
        if key == "U":
            return f["u"] + f["v"]   # the writer puts U in u (horizontal pores) or v (others)
        return f[key]

    def x(self, pore):
        L = self.net["pores"][pore]["length"]
        n = len(self.frames[0][pore]["C_bar"])
        return (np.arange(n) + 0.5) / n * L

    def path(self, pores, name, k):
        """Concatenate a field along a chain of pores, flipping pores walked outlet->inlet.
        Returns x, values, joint positions, node sequence."""
        ends = [self.net["pores"][p]["ends"] for p in pores]
        if len(pores) == 1:
            seq = [ends[0][0], ends[0][1]]
        else:
            shared = set(ends[0]) & set(ends[1])
            if not shared:
                raise ValueError(f"pores {pores[0]} and {pores[1]} do not share a node")
            start = ends[0][0] if ends[0][1] in shared else ends[0][1]
            seq = [start]
            for p, e in zip(pores, ends):
                if seq[-1] not in e:
                    raise ValueError(f"pore {p} (ends {e}) does not continue the path at node {seq[-1]}")
                seq.append(e[1] if e[0] == seq[-1] else e[0])
        xs, vs, joints, x0 = [], [], [], 0.0
        for p, e, a in zip(pores, ends, seq[:-1]):
            x = self.x(p)
            v = self.field(k, p, name)
            L = self.net["pores"][p]["length"]
            if e[0] != a:            # walked backwards
                x, v = L - x[::-1], v[::-1]
                if FIELD_ALIASES.get(name) == "U":
                    v = -v
            xs.append(x0 + x); vs.append(v)
            x0 += L; joints.append(x0)
        return np.concatenate(xs), np.concatenate(vs), joints[:-1], seq
