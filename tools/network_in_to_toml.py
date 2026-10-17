#!/usr/bin/env python3
"""Convert a legacy Network.in into network.toml (standard library only).

Usage:
    python3 tools/network_in_to_toml.py input_files/Network.in -p 1 -r 2 [-o input_files/network.toml]

Network.in has no header, so the pore and reservoir counts must be given (they are
the same numbers as num_pores / num_reservoirs in the run config). The file is read as
a whitespace token stream exactly like Network::Read_filename_1 (operator>>), so
line breaks don't matter and anything after the connectivity rows (e.g. geometry
lines for the Tecplot/Blender writer) is ignored, as the solver ignores it.

Legacy layout and how it maps:
  pore rows (10):  type dead_end in_res out_res length n_cells sigma* lambda_D area C0
  reservoir rows (9): type P_flag Mu_flag P Mu C0 sigma* lambda_D volume
  connectivity rows (num_pores): i j pore  -> dropped; TOML derives it from pore 'ends'

Old (v0) format, from the private phd_codes solver: use --v0 [--channels Channels.in]
  line 1: artificial-diffusion switch (ignored; the current solver has none)
  pore rows (4):      in_res out_res length n_cells
  reservoir rows (8): P_flag Mu_flag P Mu C0 sigma* lambda_D volume
  connectivity rows:  i j pore
  geometry rows (num_pores + num_reservoirs, 10 values): Tecplot/Blender writer, ignored
  Channels.in: per cell "sigma* lambda_D Sp C0", read as one stream in pore order,
  exactly like Read_filename_2 (rows beyond the total cell count are ignored).
  Mapping: pore type = slit (the old solver only used the slit tables), dead_end = false;
  a pore whose cells are all identical gets plain sigma_star/lambda_D/area/c0 fields,
  otherwise a per-pore CSV profile next to the output file.
  The old solver set reservoir mu = Mu + ln(C0) and P = P - 2*C0; the current solver
  uses Mu as given and forces P = 0. The converter writes potential = Mu + ln(C0) (same
  value as the old run); the pressure offset cannot be reproduced (ignored by the solver).

Quirks preserved on purpose (behavior must not change):
  * The P/Mu flags are 0 = KNOWN (fixed value), 1 = solved for, i.e. the opposite of
    input_files/README.md. The TOML uses pressure_known = (flag == 0).
  * The reservoir pressure value is read but overwritten to 0 in the solver. It is
    written to TOML (so nothing is lost) with a comment saying it is ignored.
  * The solver builds connectivity from the explicit rows; TOML uses ends=[i, j] with
    pore k. The converter checks they agree and stops if they don't.
"""
import argparse
import math
import sys
from collections import Counter
from pathlib import Path

INT_KEYS = {"n_cells"}
TYPES = {0: "slit", 1: "circular"}


def num(tok, what):
    try:
        return float(tok)
    except ValueError:
        sys.exit(f"error: expected a number for {what}, got {tok!r}")


def legacy_int(tok, what):
    # operator>> into int stops at a '.', so "1.00" would break the C++ reader too
    try:
        return int(tok)
    except ValueError:
        sys.exit(f"error: expected an integer for {what}, got {tok!r}")


def fmt(x):
    """Shortest TOML float that round-trips to the same double."""
    s = repr(float(x))
    if s in ("inf", "-inf", "nan"):
        sys.exit("error: inf/nan in Network.in")
    return s


def read_legacy(path, n_pores, n_res):
    toks = Path(path).read_text().split()
    pos = 0

    def take(n, what):
        nonlocal pos
        if pos + n > len(toks):
            sys.exit(f"error: {path} ended early while reading {what}")
        out = toks[pos:pos + n]
        pos += n
        return out

    pores = []
    for i in range(n_pores):
        t = take(10, f"pore row {i}")
        pores.append(dict(
            type=legacy_int(t[0], "pore type"), dead_end=legacy_int(t[1], "dead_end"),
            ends=[legacy_int(t[2], "inlet"), legacy_int(t[3], "outlet")],
            length=num(t[4], "length"), n_cells=legacy_int(t[5], "n_cells"),
            sigma_star=num(t[6], "sigma*"), lambda_D=num(t[7], "lambda_D"),
            area=num(t[8], "area"), c0=num(t[9], "C0")))
    res = []
    for i in range(n_res):
        t = take(9, f"reservoir row {i}")
        res.append(dict(
            type=legacy_int(t[0], "reservoir type"),
            p_flag=legacy_int(t[1], "P flag"), mu_flag=legacy_int(t[2], "Mu flag"),
            pressure=num(t[3], "P"), potential=num(t[4], "Mu"), c0=num(t[5], "C0"),
            sigma_star=num(t[6], "sigma*"), lambda_D=num(t[7], "lambda_D"),
            volume=num(t[8], "volume")))
    conn = [tuple(legacy_int(x, "connectivity") for x in take(3, f"connectivity row {i}"))
            for i in range(n_pores)]
    extra = len(toks) - pos
    return pores, res, conn, extra


def read_v0(path, n_pores, n_res, channels=None, out_dir=None, stem="pore"):
    toks = Path(path).read_text().split()
    pos = 0

    def take(n, what):
        nonlocal pos
        if pos + n > len(toks):
            sys.exit(f"error: {path} ended early while reading {what} (v0 format)")
        out = toks[pos:pos + n]
        pos += n
        return out

    take(1, "artificial-diffusion switch")
    pores = []
    for i in range(n_pores):
        t = take(4, f"pore row {i}")
        pores.append(dict(type=0, dead_end=0, ends=[legacy_int(t[0], "inlet"), legacy_int(t[1], "outlet")],
                          length=num(t[2], "length"), n_cells=legacy_int(t[3], "n_cells"),
                          sigma_star=None, lambda_D=None, area=1.0, c0=1.0))
    res = []
    for i in range(n_res):
        t = take(8, f"reservoir row {i}")
        c0 = num(t[4], "C0")
        res.append(dict(type=0, p_flag=legacy_int(t[0], "P flag"), mu_flag=legacy_int(t[1], "Mu flag"),
                        pressure=num(t[2], "P") - 2 * c0, potential=num(t[3], "Mu") + math.log(c0),
                        c0=c0, sigma_star=num(t[5], "sigma*"), lambda_D=num(t[6], "lambda_D"),
                        volume=num(t[7], "volume")))
    conn = [tuple(legacy_int(x, "connectivity") for x in take(3, f"connectivity row {i}"))
            for i in range(n_pores)]
    notes = [f"ignored {len(toks) - pos} trailing tokens (geometry rows)"] if len(toks) > pos else []
    if channels is None:
        sys.exit("error: v0 networks need --channels Channels.in (per-cell sigma*, lambda_D, Sp, C0)")
    ct = Path(channels).read_text().split()
    need = 4 * sum(p["n_cells"] for p in pores)
    if len(ct) < need:
        sys.exit(f"error: {channels} has {len(ct) // 4} cell rows, need {need // 4}")
    if len(ct) > need:
        notes.append(f"{channels}: ignored {(len(ct) - need) // 4} extra cell rows (old reader ignored them too)")
    k = 0
    for i, p in enumerate(pores):
        rows = []
        for _ in range(p["n_cells"]):
            rows.append([num(x, "Channels.in value") for x in ct[k:k + 4]])
            k += 4
        if all(r == rows[0] for r in rows):
            p["sigma_star"], p["lambda_D"], p["area"], p["c0"] = rows[0]
        else:
            name = f"{stem}{i}_cells.csv"
            with open(Path(out_dir) / name, "w") as f:
                f.write("sigma_star,lambda_D,Sp,C0\n")
                for r in rows:
                    f.write(",".join(fmt(v) for v in r) + "\n")
            p["profile"] = name
            # pore-level values are still required by the loader; use the first cell
            p["sigma_star"], p["lambda_D"], p["area"], p["c0"] = rows[0]
            notes.append(f"pore {i}: non-uniform cells -> {name}")
    return pores, res, conn, notes


def check(pores, res, conn):
    derived = sorted((p["ends"][0], p["ends"][1], k) for k, p in enumerate(pores))
    if sorted(conn) != derived:
        sys.exit("error: connectivity rows disagree with pore inlet/outlet columns:\n"
                 f"  rows:    {sorted(conn)}\n  derived: {derived}\n"
                 "fix Network.in first (TOML keeps a single source of truth: pore ends)")
    for k, p in enumerate(pores):
        if p["type"] not in TYPES:
            sys.exit(f"error: pore {k} type {p['type']} is not 0/1")
    for k, r in enumerate(res):
        if r["type"] not in TYPES:
            sys.exit(f"error: reservoir {k} type {r['type']} is not 0/1")


# Built-in defaults of the C++ loader (network_input.cpp). Precedence there is
# entry field > [defaults] > built-in, so the converter omits anything that would
# resolve to the same value anyway.
PORE_BUILTIN = {"type": 0, "dead_end": 0, "area": 1.0, "c0": 1.0}
RES_BUILTIN = {"type": 0, "pressure": 0.0, "c0": 1.0, "volume": 0.0}


def to_toml(pores, res, src):
    keys = ["type", "dead_end", "length", "n_cells", "sigma_star", "lambda_D", "area", "c0"]
    # [defaults] = values shared by every pore (skipped if equal to a built-in);
    # with a single pore everything stays on the pore itself.
    defaults = {}
    if len(pores) > 1:
        defaults = {k: pores[0][k] for k in keys
                    if all(p[k] == pores[0][k] for p in pores) and PORE_BUILTIN.get(k) != pores[0][k]}

    def effective(k):
        return defaults[k] if k in defaults else PORE_BUILTIN.get(k)

    def val(k, v):
        if k == "type":
            return f'"{TYPES[v]}"'
        if k == "dead_end":
            return "true" if v else "false"
        if k in INT_KEYS:
            return str(v)
        return fmt(v)

    L = [f"# network.toml converted from {src} by tools/network_in_to_toml.py",
         "# Schema and built-in defaults: input_files/network_schema.toml. All quantities nondimensional.",
         "# Omitted keys use [defaults], then built-ins (pores: slit, not dead-end, area 1, c0 1;",
         "# reservoirs: slit, c0 1, volume 0, pressure 0).",
         ""]
    if defaults:
        L.append("[defaults]   # used by any pore that doesn't set the key itself")
        for k in keys:
            if k in defaults:
                L.append(f"{k} = {val(k, defaults[k])}")
        L.append("")
    for k, p in enumerate(pores):
        L.append("[[pores]]")
        L.append(f"id = {k}")
        L.append(f"ends = [{p['ends'][0]}, {p['ends'][1]}]   # [inlet reservoir, outlet reservoir]")
        if p.get("profile"):
            L.append(f'profile = "{p["profile"]}"   # per-cell values override the pore values below')
        for key in keys:
            if p[key] != effective(key):
                L.append(f"{key} = {val(key, p[key])}")
        L.append("")
    for k, r in enumerate(res):
        L.append("[[reservoirs]]")
        L.append(f"id = {k}")
        if r["type"] != RES_BUILTIN["type"]:
            L.append(f'type = "{TYPES[r["type"]]}"')
        L.append(f"pressure_known = {'true' if r['p_flag'] == 0 else 'false'}")
        L.append(f"potential_known = {'true' if r['mu_flag'] == 0 else 'false'}")
        if r["pressure"] != RES_BUILTIN["pressure"]:
            L.append(f"pressure = {fmt(r['pressure'])}   # ignored: solver sets reservoir pressure to 0")
        L.append(f"potential = {fmt(r['potential'])}")
        if r["c0"] != RES_BUILTIN["c0"]:
            L.append(f"c0 = {fmt(r['c0'])}")
        for key in ("sigma_star", "lambda_D"):   # inherited from [defaults] if equal
            if defaults.get(key) != r[key]:
                L.append(f"{key} = {fmt(r[key])}")
        if r["volume"] != RES_BUILTIN["volume"]:
            L.append(f"volume = {fmt(r['volume'])}")
        L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("network_in")
    ap.add_argument("-p", "--num-pores", type=int, required=True)
    ap.add_argument("-r", "--num-reservoirs", type=int, required=True)
    ap.add_argument("-o", "--output", help="default: same folder, network.toml")
    ap.add_argument("--v0", action="store_true", help="old phd_codes format (see above)")
    ap.add_argument("--channels", help="v0 only: Channels.in with per-cell values")
    a = ap.parse_args()
    out = Path(a.output) if a.output else Path(a.network_in).with_name("network.toml")
    if a.v0:
        pores, res, conn, notes = read_v0(a.network_in, a.num_pores, a.num_reservoirs,
                                          a.channels, out.parent, out.stem + "_pore")
    else:
        pores, res, conn, extra = read_legacy(a.network_in, a.num_pores, a.num_reservoirs)
        notes = [f"ignored {extra} trailing tokens after connectivity (the solver ignores them too)"] if extra else []
    check(pores, res, conn)
    out.write_text(to_toml(pores, res, Path(a.network_in).name + (" (v0 format)" if a.v0 else "")))
    print(f"wrote {out} ({len(pores)} pores, {len(res)} reservoirs)")
    for n in notes:
        print("  note:", n)


if __name__ == "__main__":
    main()
