#!/usr/bin/env python3
"""Scale the driving potentials of the phd_codes H-junction (400 and 320) by factors and run the unchanged
solver on a given mesh; prints one JSON line per run (see tools/mesh_study.py for the fields).

    python3 tools/voltage_sweep.py --work /tmp/vsweep --mesh f1|orig|f8 --factors 1,0.5,0.25 [--dt 1e-5] [--tmax 0.03]
"""
import argparse, json, os, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mesh_study as ms  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/vsweep"); ap.add_argument("--mesh", default="orig")
    ap.add_argument("--factors", required=True); ap.add_argument("--dt", type=float, default=1e-5)
    ap.add_argument("--tmax", type=float, default=0.03); ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    cells = {0: 200, 1: 200, 2: 10, 3: 100, 4: 100} if a.mesh == "orig" else \
        {k: int(v) for k, v in ms.hj(float(a.mesh[1:])).items()}

    def one(f):
        name = f"hjV{f:g}"
        ms.CASES[name] = dict(src="examples/hjunction_phd_codes/network.toml",
                              sub=None, np=5, nr=6, nh=5, tmax=a.tmax)
        src = open(os.path.join(ms.ROOT, ms.CASES[name]["src"])).read()
        txt = src.replace("potential = 400.0", f"potential = {400*f:.6g}").replace(
            "potential = 320.0", f"potential = {320*f:.6g}")
        p = os.path.join(a.work, f"{name}_{a.mesh}.toml"); os.makedirs(a.work, exist_ok=True)
        open(p, "w").write(txt)
        ms.CASES[name]["src"] = p
        row = ms.run(a.work, name, a.mesh, cells, a.dt)
        row["factor"] = f; row["V_L"] = 400 * f; row["V_R"] = 320 * f
        return row
    with ThreadPoolExecutor(a.jobs) as ex:
        for row in ex.map(one, [float(x) for x in a.factors.split(",")]):
            print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
