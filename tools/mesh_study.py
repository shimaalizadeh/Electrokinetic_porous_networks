#!/usr/bin/env python3
"""Mesh-refinement study with the unchanged (central) solver: finds the coarsest stable
n_cells per pore for the reference H-junction, the T-junction at V=80 and the single pore.

    python3 tools/mesh_study.py --work /tmp/mesh CASE:LABEL ...   (see CASES below)
Stable = runs to max_time without "Modeling Failed" and final I_tot within 0.5% of the finest mesh.
Needs build/ek_solver and the ekpn_viz package. Appends one JSON line per run to WORK/results.jsonl.
"""
import argparse, json, os, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ekpn_viz.io import read_log, read_snapshots  # noqa: E402

EXE = os.path.join(ROOT, "build", "ek_solver")


def set_cells(text, cells):
    """Return network.toml text with n_cells set per pore id (dict id -> n)."""
    out, blocks = [], re.split(r"(?m)^(?=\[)", text)
    for b in blocks:
        if b.startswith("[[pores]]"):
            pid = int(re.search(r"(?m)^id\s*=\s*(\d+)", b).group(1))
            if pid in cells:
                b = re.sub(r"(?m)^n_cells\s*=.*\n", "", b)
                b = re.sub(r"(?m)^(id\s*=\s*\d+.*\n)", rf"\g<1>n_cells = {cells[pid]}\n", b, count=1)
        out.append(b)
    return "".join(out)


def hj(f, bridge=10):
    # dx = 1/(200 f) in the four micro-channels, bridge (pore 2) gets `bridge` cells
    n03 = round(1.3241 * 200 * f)
    return {0: 200 * f, 1: 200 * f, 2: bridge, 3: n03, 4: n03}


CASES = {
    "hj": dict(src="examples/hjunction_reference/network.toml", sub=None, np=5, nr=6, nh=5, tmax=0.03),
    "tj": dict(src="examples/tjunction/network.toml", sub=("potential = 20.0", "potential = 80.0"),
               np=3, nr=4, nh=3, tmax=0.05),
    "sp": dict(src="examples/single_pore/network.toml", sub=None, np=1, nr=2, nh=1, tmax=0.1),
}


def run(work, case, label, cells, dt, nsnap=200):
    c = CASES[case]
    d = os.path.join(work, case, label)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(os.path.join(d, "input_files"))
    os.symlink(os.path.join(ROOT, "tables"), os.path.join(d, "tables"))
    txt = open(os.path.join(ROOT, c["src"])).read()
    if c["sub"]:
        txt = txt.replace(*c["sub"])
    txt = set_cells(txt, cells)
    open(os.path.join(d, "input_files", "network.toml"), "w").write(txt)
    we = max(1, round(c["tmax"] / dt / nsnap))
    open(os.path.join(d, "run.toml"), "w").write(
        f'[network]\nnum_pores = {c["np"]}\nnum_reservoirs = {c["nr"]}\nnum_horizontal_pores = {c["nh"]}\n'
        f'[solver]\ndt = {dt}\nmax_time = {c["tmax"]}\nwrite_every = {we}\n'
        f'[io]\ntables_dir = "tables/"\nnetwork_file = "input_files/network.toml"\n')
    t0 = time.time()
    with open(os.path.join(d, "log.txt"), "w") as fh:
        rc = subprocess.run([EXE, "run.toml"], cwd=d, stdout=fh, stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    log = open(os.path.join(d, "log.txt")).read()
    t, I = read_log(os.path.join(d, "log.txt"))
    row = dict(case=case, label=label, cells=cells, dt=dt, wall_s=round(wall, 1),
               failed=bool("Modeling Failed" in log or rc != 0 or len(t) == 0 or t[-1] < c["tmax"] * 0.999),
               t_end=float(t[-1]) if len(t) else 0.0, I_end=float(I[-1]) if len(I) else None)
    try:
        _, frames = read_snapshots(d, dt, we)
        row["minC"] = float(min(min(fr[p]["C_bar"].min() for p in range(len(fr))) for fr in frames))
        if case == "hj":
            row["bridge_minC"] = float(min(fr[2]["C_bar"].min() for fr in frames))
    except Exception as e:  # noqa: BLE001
        row["minC_err"] = str(e)[:80]
    for fn in os.listdir(d):  # keep disk use small
        if fn.startswith("Tecplot_"):
            os.remove(os.path.join(d, fn))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/mesh")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("specs", nargs="+", help="case:label:dt:cells, cells = f<F>[b<B>] for hj, N (uniform) or n0,n1,... (per pore)")
    a = ap.parse_args()
    jobs = []
    for s in a.specs:
        case, label, dt, cs = s.split(":")
        if "," in cs:  # explicit per-pore list n0,n1,...
            cells = {i: int(v) for i, v in enumerate(cs.split(","))}
        elif case == "hj":
            m = re.fullmatch(r"f([\d.]+)(?:b(\d+))?", cs)
            cells = hj(float(m.group(1)), int(m.group(2) or 10))
            cells = {k: int(v) for k, v in cells.items()}
        else:
            n = int(cs); cells = {p: n for p in range(CASES[case]["np"])}
        jobs.append((case, label, cells, float(dt)))
    with ThreadPoolExecutor(a.jobs) as ex:
        for row in ex.map(lambda j: run(a.work, *j), jobs):
            print(json.dumps(row), flush=True)
            with open(os.path.join(a.work, "results.jsonl"), "a") as fh:
                fh.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
