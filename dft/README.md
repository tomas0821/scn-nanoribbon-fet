# DFT and Wannier90 data

All DFT runs used Quantum ESPRESSO 7.5 and Wannier90. Every stage directory keeps the files that generate, run and
parse its calculations, plus the inputs and the parsed results:

| Kind | Files |
|---|---|
| input generators | `make_*.py` (write `inputs/*.in`, the Wannier90 `.win` files and `tags*.txt`, the ordered run lists) |
| parsers | `parse_*.py` (read QE / Wannier90 outputs and write the stage's `*.json` summaries) |
| job scripts | `run_*.sh` (SLURM scripts as used; the site-specific partition, account and QOS lines are generalised) |
| QE inputs | `inputs/*.in` (`pw.x`, `projwfc.x`, `pp.x`, `open_grid.x`) |
| Wannier90 inputs and logs | `w90/**/<seed>.win`, `<seed>.pw2wan`, `<seed>.wout` (final centres and spreads), `<seed>_centres.xyz` |
| Wannier Hamiltonians | `w90/**/<seed>_hr.dat` (only the ones the analysis scripts read; see below) |
| orbital maps | `w90/<ribbon>/wf_map.json`: atom, species and orbital of every Wannier function, per projection set (`B`: Sc d + N p (+ H s, + Sc s on the [110] bare edge Sc); `S`: Sc s on every Sc) |
| summaries | `*.json`: band edges, vacuum levels, Wannier windows, pass/fail criteria, costs |
| raw outputs | `outputs/*.out`: only those the analysis scripts parse (relaxed geometries, band and grid eigenvalues) |

Paths to the pseudopotential directory and to the cluster scratch area read `/path/to/scn_fet/...`; set them to a
local directory before rerunning. The pseudopotentials are the published files named in the inputs:
`Sc.pbe-spn-kjpaw_psl.0.2.3.UPF` (SSSP precision, PAW) and `N.oncvpsp.upf` (PseudoDojo) for the PBE relaxations;
`Sc_ONCV_PBE-1.0.oncvpsp.upf` and `H_ONCV_PBE-1.0.oncvpsp.upf` (SG15) with `N.oncvpsp.upf` for the norm-conserving
PBE-twin and HSE06 runs.

## Stages

Each stage directory also has its own `README.md`, kept as written during the project apart from wording (decision
and approval notes as a neutral, dated record; pointers to undeposited project files removed): the registered criteria
and their dated amendments, the run log (job numbers, core-hours) and the results. Lines 1–84 of `A3/README.md` are
checksummed in `A3/README.prereg.md5`, which gives the deposited copy's hash first and the frozen original's after it
(see "Pre-registration" in the top-level README). The labels used there (A1-1, Gate M, P1–P6, H1–H4, LR3, …) are
listed in the top-level README. Those READMEs give cluster paths as `/scratch/global/<user>/scn_fet`; the input
files here read `/path/to/scn_fet` instead.

| Stage | What | Key files |
|---|---|---|
| `A0` | cutoff and pseudopotential tests (bulk ScN and the M-ScN sheet) | `a0_summary.json` |
| `A1` | bulk rocksalt ScN: PBE and HSE06 bands, Wannier90 on a 16³ grid (`w90/hse_k16`, `w90/pbe_k16`; the other seeds are earlier diagnostics) | `w90/{hse,pbe}_k16/*_hr.dat`, `outputs/{hse_ibz,pbe_nscf}_k16.out` (SK fits) |
| `A2` | the M-ScN sheet (two-plane rocksalt-(001) bilayer): PBE and HSE06, Wannier90 on a 16 × 16 grid; `costtest/` holds the HSE cost and band-count tests | `w90/{hse,pbe}_k16/` (building-block interior, SK fit), `a2_stage1.json`, `a2_w90.json` (vacuum levels) |
| `A3` | three-plane slab hold-out: PBE relaxation, NC PBE twin + Wannier90, HSE06 + Wannier90 | `outputs/relax_slab3.out`, `w90/slab3/W2`, `w90/slab3/H/W2`, `a3_criteria.json`, `a3_hwindows.json` |
| `A4a` | calibration ribbons `r100_3`, `r110_4`, `r100_3h`, `r110_4h` ([100] n = 3 and [110] n = 4, bare and H-terminated): PBE relaxations, magnetism and vacuum checks | `outputs/relax_<ribbon>.out` |
| `A4b` | calibration ribbons: NC PBE twin + Wannier90 (A4b-1) and HSE06 + Wannier90 (A4b-2, under `w90/<ribbon>/H/`) | `w90/<ribbon>/{W1,W2}`, `w90/<ribbon>/H/{W1,W2,W3}`, `a4b_windows.json`, `a4b2_windows.json`, `outputs/{bands,hse}_<ribbon>.out` |
| `A5` | hold-out ribbons `r110_8` and `r100_6`: PBE relaxations (`outputs/`), NC PBE twin + Wannier90 (`twin/`), HSE06 gap single points (`hse/`) | `outputs/relax_<ribbon>.out`, `twin/w90/<ribbon>/W2`, `hse/a5_hse.json` |

Settings, in brief: PBE relaxations with the SSSP-precision PAW set at 90/720 Ry; norm-conserving PBE twins (smearing:
Gaussian for the ribbons, Marzari–Vanderbilt for the slab) and HSE06 (fixed occupations, `ecutfock` 160 Ry) at 80 Ry; ribbons periodic along x with 15 Å ([100] bare) or 20 Å
([110] and H-terminated) of vacuum. Apart from the generalised paths, each input file is the exact input that was run.

## Wannier variants

Each ribbon and the slab have several Wannier90 runs that differ only in the disentanglement windows and projection
set: `W1` (the sheet's window rule), `W2` (`W1` with the frozen window stopped 0.1 eV below the first vacuum-like
state), `W3` (projection set `S`, an Sc s function on every Sc). The prefix `H` marks the HSE06 runs (`HW1`, `HW2`,
`HW3`, stored under `H/`). The device Hamiltonians of the calibration ribbons are `r100_3` HW3, `r110_4` HW1,
`r100_3h` HW2 and `r110_4h` HW2; the building blocks use the edges of W1 / HW1. The `.win` and `.wout` files of every
variant are included; the `_hr.dat` files of variants that no script reads are not.
