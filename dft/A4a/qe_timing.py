"""Timing facts from our own QE 7.5 outputs (A0–A2), for the A4a/A4b/A5 cost model (cost_a4.py)."""
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _f(pat, txt, cast=float, default=None):
    m = re.search(pat, txt)
    return cast(m.group(1)) if m else default


def wall_s(txt, routine):
    """Wall seconds of a timing line, e.g. 'vexx : 3709.66s CPU 3809.43s WALL ( 36 calls)' (h/m forms too)."""
    m = re.search(rf"^\s*{routine}\s*:\s*(.*?)CPU\s+(.*?)WALL(?:\s*\(\s*(\d+) calls\))?", txt, re.M)
    if not m:
        return None, None
    s = m.group(2).strip()
    t = 0.0
    for v, u in re.findall(r"([\d.]+)([hms])", s):
        t += float(v) * dict(h=3600, m=60, s=1)[u]
    return t, (int(m.group(3)) if m.group(3) else None)


def parse(out, inp=None):
    txt = Path(out).read_text(errors="replace")
    d = dict(name=Path(out).stem, stage="/".join(Path(out).parts[-4:-2]).removeprefix("dft/"))
    d["ranks"] = _f(r"Number of MPI processes:\s+(\d+)", txt, int, 1)
    d["npool"] = _f(r"npool\s*=\s*(\d+)", txt, int, 1)
    d["nbnd"] = _f(r"number of Kohn-Sham states=\s*(\d+)", txt, int)
    d["nel"] = _f(r"number of electrons\s*=\s*([\d.]+)", txt)
    d["nks"] = _f(r"number of k points=\s*(\d+)", txt, int)
    d["nspin"] = 2 if "LSDA" in txt or "spin-polarized" in txt.lower() else 1
    d["ecut"] = _f(r"kinetic-energy cutoff\s*=\s*([\d.]+)", txt)
    d["vol_au"] = _f(r"unit-cell volume\s*=\s*([\d.]+)", txt)
    d["nat"] = _f(r"number of atoms/cell\s*=\s*(\d+)", txt, int)
    m = re.search(r"Sum\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)", txt)
    d["npw_sum"] = int(m.group(1)) if m else None              # smooth-grid PWs (total over G)
    d["iters"] = len(re.findall(r"iteration #\s*\d+", txt))
    d["wall"], _ = wall_s(txt, "PWSCF")
    vx, calls = wall_s(txt, "vexx")
    if vx:
        g = re.search(r"EXX grid:\s+\d+ G-vectors\s+FFT dimensions: \(\s*(\d+),\s*(\d+),\s*(\d+)\)", txt)
        d["Nr"] = int(g.group(1)) * int(g.group(2)) * int(g.group(3))
        d["vexx"], d["vexx_calls"] = vx, calls
        d["outer"] = len(re.findall(r"EXX: now go back to refine", txt)) + 1
        if inp:
            it = Path(inp).read_text()
            nq = [int(_f(rf"nqx{i}\s*=\s*(\d+)", it, int, 1)) for i in (1, 2, 3)]
            d["nq"] = nq[0] * nq[1] * nq[2]
    return d


def hse_rows():
    rows = []
    for st in ("A1", "A2", "A2/costtest"):
        for out in sorted((ROOT / st / "outputs").glob("hse_*.out")) + sorted((ROOT / st / "outputs").glob("ct_*.out")):
            if out.stem.endswith("_open") or "pot" in out.stem:
                continue
            inp = ROOT / st / "inputs" / f"{out.stem}.in"
            d = parse(out, inp if inp.exists() else None)
            if d.get("vexx") and d.get("nq"):
                rpp = d["ranks"] // d["npool"]
                kpp = math.ceil(d["nks"] / d["npool"])                 # rank 0's pool holds the most k
                nocc = int(round(d["nel"])) // 2
                pairs = d["vexx_calls"] * d["nq"] * d["nbnd"] * nocc   # per process (one pool)
                d["C_pair"] = d["vexx"] * rpp / pairs                   # core-s per (k, q, i, j) pair
                d["C_norm"] = d["C_pair"] / (d["Nr"] * math.log(d["Nr"]))
                d["C_sq"] = d["vexx"] * rpp / (d["vexx_calls"] * d["nq"] * d["nbnd"] ** 2 * d["Nr"] * math.log(d["Nr"]))
                d["outer_from_calls"] = d["vexx_calls"] / kpp
                rows.append(d)
    return rows


if __name__ == "__main__":
    for d in hse_rows():
        print(f"{d['stage']} {d['name']:14s} nks {d['nks']:3d} npool {d['npool']} nq {d['nq']:4d} nbnd {d['nbnd']:3d} "
              f"nocc {int(d['nel'])//2:2d} Nr {d['Nr']:7d} outer {d['outer']} (calls/kpp {d['outer_from_calls']:.1f}) "
              f"vexx {d['vexx']:7.0f}s wall {d['wall']:7.0f}s  C_pair {1e3*d['C_pair']:.3f} ms  "
              f"C_norm {1e9*d['C_norm']:.3f} ns  C_sq {1e9*d['C_sq']:.2f} ns")
