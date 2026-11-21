# A5: the two pre-registered hold-out ribbons (plan §4 A5, §3 D7)

**Scope (decision 2026-11-14, option 2):** the plan's original A5.
- PBE relaxation, PBE twin with Wannier90 (the "PBE-Wannier hold-outs" of D7), and HSE06 **gap single points**.
- No HSE-Wannier device Hamiltonians for the hold-outs. That reverses LR1 option A for these two ribbons only, after
  the A4b-2 recalibration put HSE-W90 on them at 4.8–12.7k core-h with 34–108 h per one-node run.
- The four calibration ribbons keep their HSE device Hamiltonians (A4b-2).
- The hold-outs' device results come from the SK model, tested here as D7 says.

**Written 2026-11-14, before any A5 input or output exists.** Later changes are dated amendments below, never edits.

## Ribbons (never refitted)
| Key | Ribbon | Atoms | W (edge to edge) | Vacuum | Cell (Å) |
|---|---|---|---|---|---|
| `r110_8` | [110] n = 8, bare | 32 (16 Sc) | 2.33 nm | 20 Å (decision 12 Nov: [110] at 20 Å) | 3.109 × 43.32 × 22.17 |
| `r100_6` | [100] n = 6, bare | 48 (24 Sc) | 2.42 nm | 15 Å (decision 13 Nov: bare [100] 0.0 meV at 15 → 20 Å) | 4.397 × 39.19 × 17.17 |

Built by `scnfet.structures.ribbon` with the A2/A4a convention (a_p 3.1093 Å, planes 2.1685 Å, fixed period).

## Stages
**A5-1, PAW PBE, as A4a** (Sc PAW psl 0.2.3 spn, N Dojo; 90/720 Ry; gaussian 0.005 Ry; k 15):
1. `relax_<r>`: internal coordinates from the ideal cut, nspin 1, BFGS, `forc_conv_thr` 1e-3, `etot_conv_thr` 1e-5.
   - Non-magnetic relaxation because Gate M found both edge types non-magnetic (every magnetic start collapsed at
     n = 3 and 4).
   - The magnetism check follows in item 2.
   - One job per ribbon. `max_seconds` + `restart_mode='restart'` as in A4a.
2. At the relaxed geometry (staged after item 1):
   - `nm_<r>`, `fm_<r>`, `afmedges_<r>` single points at k 15;
   - `bands_<r>` on the 31-k Γ–X̄ path, with `projwfc_<r>` on it.
   - **Magnetism, Gate M criterion 1 as registered:** magnetic if ΔE_mag > 5 meV per edge atom. If magnetic →
     plan LR5, back to the authors.
   - The band-edge character (edge weights, as `edge_bands.py`) is descriptive.

**A5-2, NC PBE twin, as A4b-1** (`dft/A4b/README.md`, recipe frozen 13 Nov; set B with Sc s on the edge Sc of bare
[110]):
- scf, pot, bands, projwfc, NSCF on the Wannier grid (16 k [110] / 12 k [100]), `pw2wannier90`, then W1/W2 (and W3
  only where amendment 1's rule applies).
- **Criteria:**
  - P1–P6 as frozen;
  - the decision rule W2 → W1 → back to the authors (LR6);
  - **plus H4** (amendment 2, pre-registered for the hold-outs): every pair of variants agrees within 5 meV over the
    device window on a 601-point path, otherwise back to the authors.
- **To expect:** in r110_4 the CB edge pair scored s_out 0.248 in PBE, just under 0.25. If it crosses 0.25 in
  r110_8, W2's frozen top falls below the CB and the rule falls back to W1. That is handled by the rule, not tuned.
- The Wannier90 Hamiltonians are the PBE-Wannier hold-outs that the PBE-twin SK model is tested against (D7).

**A5-3, HSE06 gap single points, as A4b-2** (NC, 80 Ry, `ecutfock` 160, fixed occupations; IBZ of the 16 / 12-k grid;
q 8 / 6):
- **nbnd = nocc + 32** (160 for r110_8, 224 for r100_6). The gap needs the band edges only; no Wannier outer window.
- **Also run:** `projwfc.x` on the HSE run and `pp.x` for the vacuum level.
- **The vacuum-like state is reported, not a stop:** no HSE device Hamiltonian is built here. If it lies inside
  the device window (VBM − 0.5 to CBM + 0.5 eV), the SK device results at that width carry it as a stated limitation,
  because the SK basis (Sc d + N p) has no such band.
- **Gap:** from the IBZ grid eigenvalues, with the `!!` convergence check of `parse_a4b.py`.

## Pass criteria (plan §3 D7, pre-registered in the plan; judged at the claims gate, not here)
- **PBE twin SK vs these PBE-Wannier hold-outs:**
  - gap ±50 meV;
  - valley offsets ±50 meV;
  - transport masses ±10 %;
  - I_ON ±15 %;
  - the same IRDS verdict.
- **SK-HSE gap vs the HSE gap single points:** ±80 meV.
- Failure → plan LR2: no refitting; reported as a transferability limit.

## Cost (`dft/A4a/cost_a4.py`, recalibrated 14 Nov)
| Stage | core-h | Wall |
|---|---|---|
| A5-1 relax | r110_8 213–726, r100_6 357–1,213 (law × 1.5–5.1, the A4a range) | 3–11 h / 6–19 h on one node |
| A5-1 single points, bands, projwfc | ≈ 100–300 | |
| A5-2 PBE twin + Wannier90 | ≈ 300–700 | |
| A5-3 HSE gap points (+ projwfc, pp.x) | r110_8 450–960, r100_6 480–1,030 | 7–16 h each on one node |
| **A5 total** | **≈ 1.9–4.9k** | |

The project total with A5 and A3 (≈ 100–200) comes to ≈ 7.5–10.6k core-h (spent 5,484).

## Results: A5-1 relaxations (2026-11-15)
- **Jobs 143316 (`r110_8`, cn003, 3 h 38 min) and 143317 (`r100_6`, cn004, 4 h 08 min):** both BFGS-converged in one
  job with no restart: 16 / 11 BFGS steps (17 / 12 SCF cycles), final total force 0.0022 Ry/Bohr each.
- **Cost:** 232 + 265 = **497 core-h** (estimate 213–726 and 357–1,213). `r100_6` came in below its range: the relax
  needed only 12 SCF cycles. Spent **5,981**.
- **Geometry (descriptive), against the A4a calibration ribbon with the same edge:**

  | | `r110_8` | `r110_4` (A4a) | `r100_6` | `r100_3` (A4a) |
  |---|---|---|---|---|
  | Largest displacement, edge row (Å) | 0.226 | 0.235 | 0.068 | 0.074 |
  | Largest displacement, interior (Å) | 0.150 (second row) | 0.138 | 0.038 | 0.030 |
  | Relaxation energy (eV per cell) | −2.17 | −2.05 | −0.35 | −0.34 |
  | W, ideal → relaxed (nm) | 2.332 → 2.312 | | 2.419 → 2.412 | |

  - The relaxation stays in the outer two rows. The centre rows move ≤ 0.01 Å.
  - Displacements are mirror images across the ribbon centre, so the relaxed cut keeps its symmetry.
  - Both ribbons relax like their n = 3 / 4 counterparts. The extra rows add nothing new at the edge.
- **Width:** the pre-registered W (2.33 / 2.42 nm, edge to edge between nuclei) is the ideal cut. The relaxed width is
  0.02 / 0.01 nm smaller.

## A5-1b staged (2026-11-15)
`make_a5.py stage2` → `inputs/{nm,fm,afm_edges,bands,projwfc}_<r>.in`, `tags_stage2_<r>.txt`, `a5_inputs.json`
(stage2 block).
- **What runs, as registered in item 2 above:**
  - NM, FM and AFM-across-the-edges SCFs at k 15 (A4a's `single_point` and `spin_labels`, A4a's starting moments);
  - `bands` on the 31-k Γ–X̄ path and `projwfc` on it.
- **Tag name:** the run uses `afm_edges_<r>` (A4a's spelling) for item 2's `afmedges_<r>`.
- **Not run:** A4a's AFM-between-planes start and its k 30 check. A5 did not register them; Gate M found k 15
  enough.
- **Band count of the `bands` run:** 595 for `r110_8`, 672 for `r100_6`.
  - This is A4a's rule: 28 per f.u., with the empty bands scaled by the cell's yz area beyond 15 Å (×1.46 for
    `r110_8` at 20 Å).
  - The run is also the band survey that sets A5-2's NSCF band count, as A4a's did for A4b-1.
- **Gate M criterion 1:** per edge atom, with the edge atoms counted as in Gate M (outermost row, both planes): 4 for
  `r110_8`, 8 for `r100_6`.
- **Cost: 515–744 core-h; one node, 3.2–4.7 h (`r110_8`) and 4.8–7.0 h (`r100_6`).**
  - Basis: each relaxation's first, from-scratch NM SCF, times the ratios A4a-2 measured on its ribbons: FM 1.73–1.81,
    AFM 1.67–1.70, bands 2.13–2.74, projwfc ≈ 0.08, plus 30 % at the top for SCF iteration counts. The figures are in
    `cost_a4.py` → "After the A5-1 relaxations".
  - **Correction:** this replaces the 100–300 in the cost table above, which took the cost law at face value.
  - **Net:** A5 at the relaxations' actual cost is ≈ 2.2–3.9k (was 1.9–4.9k). The project comes to ≈ 7.8–9.6k (was
    7.5–10.6k).
- **Stub dry run of `run_a5.sh` with the tag lists:** the five tags run in order, `tmp_nm_<r>` is kept for bands and
  projwfc, and the relaxations are not touched.

## Results: A5-1b (2026-11-15)
Jobs 143514 (`r110_8`, cn004, 3 h 47 min) and 143515 (`r100_6`, cn007, 5 h 26 min), every tag done. **590 core-h**
(estimate 515–744); spent **6,571**. `parse_a5.py` → `a5_stage2.json` (A4a's Gate M helpers, unchanged).

**Magnetism, Gate M criterion 1 (registered): both hold-outs are non-magnetic. No LR5.**

| | FM start | AFM-edges start | ΔE_mag (NM − lowest magnetic) |
|---|---|---|---|
| `r110_8` (4 edge atoms) | \|m\| 2.88 → 0.00 μB | 2.91 → 0.00 μB | −0.27 meV per edge atom (NM lowest) |
| `r100_6` (8 edge atoms) | \|m\| 4.37 → 0.00 μB | 4.35 → 0.00 μB | −0.24 meV per edge atom (NM lowest) |

Every magnetic start collapsed to zero moment, as at n = 3 and 4.

**PBE band edges (descriptive), on the 31-k path:**
- **`r110_8`: gap 0.264 eV, direct at Γ** (`r110_4`: 0.409).
  - The CB edge is again the Sc-s edge pair (LR3 holds at n = 8): edge weight 0.88, Sc s 0.61, |ψ|² 0.82.
  - Its two members are 2 meV apart: the edges no longer talk to each other.
  - The first interior CB state (Sc d) lies at CBM + 0.54 eV (VBM + 0.80).
  - The VB top is interior N p: edge weight 0.13, uniform share 0.125.
  - **For A5-2:** in PAW this pair would score s_out ≈ 0.23 against set B, close to the 0.25 line, as in `r110_4`
    (0.248). The NC projections of A5-2 decide it. If it crosses, the registered rule falls back to W1.
- **`r100_6`: gap 0.612 eV, direct at Γ** (`r100_3`: 0.912). Both edges are interior: CB Sc d (edge weight 0.10,
  uniform 0.167) and VB N p (edge weight 0.21).
  - A vacuum-like state (Sc s 0.31, N s 0.09, |ψ|² 0.76) sits at Γ, CBM + 0.71 eV. That is outside the device window
    in PBE.
  - It is the same kind of state as `r100_3`'s: PBE + 0.52, which HSE moved to + 0.35 and so inside the window.
  - A5-3 reports where HSE puts it (a stated SK limitation if it is inside the device window; not a stop).
- The gaps shrink with width towards the sheet's PBE 0.404 eV.

**Band survey for A5-2** (smallest N whose band minimum on the path clears the level):

| | N(VBM + 4 eV) | N(L = VBM + 10.28) | N(L + 0.5) | N(L + 1.0) |
|---|---|---|---|---|
| `r110_8` | 185 | 359 (22.4 per f.u.) | 380 | 400 |
| `r100_6` | 251 | 472 (19.7 per f.u.) | 488 | 499 |

## A5-2 staged: job 1 (2026-11-15)
`make_a5.py twin` → `twin/inputs/{scf,pot,bands,projwfc,nscf}_<r>.in`, `twin/w90/<r>/<r>_{B,S}.{win,pw2wan}`,
`twin/tags_<r>.txt`.
- **Code:** `make_a4b.ribbon_inputs`, unchanged, pointed at the A5 relaxed geometry. That geometry already sits in
  its A5 cell, and the positions match the relaxation output exactly.
- **Isolation:** the inputs live in `twin/`, so no A5-1 file is touched (its `bands_<r>` / `projwfc_<r>` names would
  collide). The 80 A4b inputs and the A5-1 inputs are unchanged (md5).
- **Settings, as A4b-1:** NC set, 80/320 Ry, gaussian 0.005 Ry, SCF at k 15, the 31-k path, NSCF on the full W90 grid
  (16 k [110], 12 k [100]) with nosym/noinv.
- **Projection sets:**
  - set B = Sc d + N p, plus Sc s on the two edge Sc of `r110_8` (as `r110_4`);
  - set S = B + Sc s on every Sc. Its `.amn` is cheap, and W3 runs only where amendment 1's rule applies.
- **Band counts, A4b-1's rules applied to the survey:**

  | | bands run | NSCF | num_bands | num_wann B / S |
  |---|---|---|---|---|
  | `r110_8` | 200 | 415 | 335 | 130 / 144 |
  | `r100_6` | 275 | 515 | 395 | 192 / 216 |

  - bands = N(VBM + 4) + 8 %; NSCF = N(L + 1) + 3 %; both rounded up to a multiple of 5.
  - A4b-1 had N(VBM + 4) + 10 bands and N(L + 1) + 8–9 bands. Excluded bands: Sc 3s 3p and N 2s (80 / 120).
- **Job:** A4b-1's `run_a4b.sh`, unchanged, run in `scn_fet/A5/twin` via `A4B_DIR`. One 64-rank job per ribbon
  (8 h limit), 8 / 4 pools on the NSCF.
- **Stub dry run:** all 8 tags per ribbon, a rerun skips all 8, and the NSCF pool counts are right. (The stub first
  read the wrong argument for `pw.x`; fixed and rerun.)
- **Cost, scaled from A4b-1's measured jobs (31 min each) by the PAW size ratios measured in A5-1b:**
  - job 1 ≈ 230–390 core-h (`r110_8` 1.5–3.1 h, `r100_6` 2.1–3.1 h);
  - job 2 (Wannier90 W1/W2, W3 only if needed) ≈ 130–190 core-h. It is staged after job 1, from
    `parse_a4b.py windows`. Its run script needs more than A4b-1's 2 h: serial Wannier90 with 130–216 WFs is ≈ 5–7×
    A4b-1's.
  - A5-2 ≈ 360–580 core-h (README estimate 300–700).

## A5-3 staged (2026-11-15)
`make_a5.py hse` → `hse/inputs/{hse,hpot,hproj}_<r>.in`, `hse/tags_hse_<r>.txt`; job script `hse/run_a5_hse.sh`.
- **Inputs:** `make_a4b.hse_inputs`, unchanged, at the relaxed geometry (the positions equal the twin's).
  - Only its first three steps are kept (HSE SCF, vacuum potential, projections). Its open_grid / Wannier90 inputs are
    removed after generation.
  - Against A4b-2's `hse_r110_4.in` the header differs only in nat and nbnd. `hpot` / `hproj` are identical to
    A4b-2's.
- **Settings as registered:** nbnd = nocc + 32 = **160** (`r110_8`) / **224** (`r100_6`); IBZ 9 / 7 k of the 16 / 12-k
  grid; q 8 / 6; fixed occupations; `ecutfock` 160.
- **Job script `run_a5_hse.sh`:** `run_a4b2.sh`'s three steps, with two changes:
  1. The rank count follows the allocation: 63 on one node, 126 on two, one pool per IBZ k either way.
  2. `hse_<r>` counts as done only with QE's converged `!!` energy. If it has none (e.g. a max_seconds stop, which
     also prints JOB DONE), the job stops before `hpot` / `hproj`.
- **Stub dry run:** one node 63 ranks / 7 pools; two nodes 126 / 9; a rerun skips all three; an unconverged HSE
  stops the job with no potential or projections written.
- **Memory:** the nodes have ≈ 1 TB. Scaling A4b-2's measured peaks (`sacct` MaxRSS 7.1–11.4 GB per process) by
  FFT slice × nbnd × k+q points gives ≈ 9–10 GB per process, ≈ 560–650 GB on one node. One node fits either ribbon.
- **Wall time:** scaling A4b-2's runs by nbnd² × EXX grid (the measured law) gives `r110_8` ≈ 11–13 h on one node.
  With the ×1.7 node spread that comes near the 22 h `max_seconds` stop, and there is no EXX restart.
  - **Proposed: `r110_8` on two nodes** (≈ 6–8 h), `r100_6` on one (≈ 7.5 h; registered range 7.5–16 h).
  - The core-h stay as registered (`r110_8` 450–960, `r100_6` 480–1,030), plus perhaps 10–20 % for two-node
    efficiency.

## Results: A5-2 job 1, the NC PBE twins (jobs 143633 / 143634, 2026-11-15)
- **Runs:** `r110_8` 1 h 37 min (cn004), `r100_6` 2 h 00 min (cn007); 231 core-h in total (estimate 230–390).
  - All 8 tags done in each; pw2wannier90 `JOB DONE` for sets B and S.
  - Copied back: `outputs/`, the pot cubes, `.nnkp` / `.eig` / `.wout` / `pw2wan.out` (not `.amn` / `.mmn`).
- **Windows: `parse_a5.py twin_windows`.** This is `parse_a4b.py windows`, unchanged, pointed at `twin/`:
  - its make_a4b is `make_a5._make_a4b(TWIN)` with r110_8's edge Sc s;
  - its run facts are in `twin/a4b_inputs.json`, written from `a5_inputs.json`.
  - Output: `twin/a4b_windows.json`, `twin/w90/<r>/{W1,W2,W3}/<r>_<v>.win`. A4b's own json files are unchanged (md5).

| | `r110_8` | `r100_6` |
|---|---|---|
| NC PBE gap (path = grid) | 0.285 eV (PAW 0.264) | 0.617 eV (PAW 0.612) |
| E_vac − VBM (spread) | 3.115 eV (23.5 meV) | 3.351 eV (5.9 meV) |
| Edge pair s_out (set B), max on the path | 0.239 (r110_4: 0.248) | — (no edge state at the edges) |
| First non-model state (s_out ≥ 0.25) | CBM + 1.105 eV, Γ, band 133 | CBM + 0.714 eV, Γ, band 196 |
| … at 0.20 / 0.30 | the edge-pair CBM / + 1.105 | + 0.714 / + 0.714 |
| W1 frozen top (bound) | CBM + 2.330 (vacuum) | CBM + 2.234 (vacuum) |
| W2 frozen top (bound) | CBM + 1.005 (non-model) | CBM + 0.614 (non-model) |
| Band counts W1 / W2 / W3 | ok / ok / ok | ok / ok / ok |
| Outer window | VBM + 10.23, not capped | VBM + 10.23, not capped |

- **The [110] edge pair stays under 0.25** (0.239), so W2 freezes above it: the registered rule holds without the
  W1 fallback.
  - It sits between the 0.20 and 0.25 lines, as in r110_4. The 0.20 sensitivity check would put W2's top at the
    CBM. That is reported, not used.
- **The first non-model state is the same in both:** an interior state, not an edge state (edge weight 0.00 /
  0.02).
  - Character: Sc d 0.60–0.62, Sc s 0.30, N s 0.08–0.09, |ψ|² 0.77.
  - This is the sheet's vacuum-reaching state (sheet PBE: CBM + 0.885, same character; `dft/A2/a2_sheet_nonmodel.json`).
  - In [100] it falls towards the CBM as the ribbon narrows: sheet + 0.885 → n = 6 + 0.714 → n = 3 + 0.522 eV.
- **Amendment 1, PBE: no non-model state in the device window** (≤ CBM + 0.5) for either ribbon, so W3 is not run
  (its `.win` files are written).
  - For [100] n = 6, the HSE run (A5-3) may move it into the window: in n = 3, HSE moved it from CBM + 0.52 to
    + 0.35.
  - If so, it is reported as the registered limitation; there is no HSE device Hamiltonian here.
- **W2 freezes the whole device window in both** (CBM + 1.005 / + 0.614 ≥ CBM + 0.5).

## A5-2 job 2 staged (2026-11-15), with the slab's PBE variants (A3-3)
`twin/run_a5_w90.sh` runs `run_a4b_w90.sh`'s steps (serial `wannier90.x` per variant, set-B `.amn/.mmn/.eig` linked).
- **Arguments are `<dir>:<ribbon>`**, so the A3 slab (`A3:slab3`) runs on the same node.
- **Variants W1 W2** (W3 only with `VARIANTS="W1 W2 W3"`, where amendment 1's rule applies: nowhere so far).
- **8 h limit.** The longest run is ≈ 1.1–3.5 h: `r100_6` with 192 WFs and 395 bands is ≈ 5–9× r100_3's 7–12 min
  per variant.
- **Stub dry run:** 6 variants started with the set-B links; a rerun skips them.
- **Cost:** one node for ≈ 1.1–3.5 h, ≈ 70–225 core-h (README: job 2 130–190).
- **After it:** `parse_a4b.py criteria` on the twin (P1–P6, decision rule) plus H4 between W1 and W2 (601-point
  path), and the 2D version for the slab.

## A5-2 job 2 submitted (2026-11-15, approved)
The script and 6 variant `.win` files were md5-checked into `scn_fet`. **Job 143673** on cn004, started 16:06:54,
running from `A5/twin`: `sbatch -J a5w run_a5_w90.sh A5/twin:r110_8 A5/twin:r100_6 A3:slab3`. After 30 s all six
variants were in the disentanglement with the outer windows as set.

## Results: A5-2 job 2, the PBE-Wannier hold-outs (job 143673, 1 h 09 min, 74 core-h; 2026-11-15)
All six Wannier90 runs ended with "All done" (both hold-outs and the slab; the slab is in `dft/A3/README.md`). Copied
back: `.wout`, `_hr.dat`, `_centres.xyz`, `_wsvec.dat`, `w90.log` (not the linked `.amn` / `.mmn`).
`parse_a5.py twin_criteria` = `parse_a4b.py criteria`, unchanged, on the twin (→ `twin/a4b_criteria.json`), then H4
between the variants (→ `twin/a5_h4.json`). A4b's own json files are unchanged (md5).

| | `r110_8` W1 | `r110_8` W2 | `r100_6` W1 | `r100_6` W2 | Limit |
|---|---|---|---|---|---|
| P1 (vacuum margin, counts, outer) | ok | ok | ok | ok | |
| P2 device-window max error | 1.27 meV | 4.71 meV | 0.19 meV | 0.72 meV | ≤ 5 |
| P3 frozen RMS (max) | 0.55 (15.8) | 1.93 (41.8) | 0.62 (11.0) | 0.08 (0.8) | ≤ 3 |
| P4 gap | −0.00 meV | −0.00 meV | −0.01 meV | +0.01 meV | ≤ 3 |
| P5 m* VB / CB | +0.2 / +1.3 % | +0.4 / +2.9 % | +0.0 / +0.3 % | −0.0 / +0.5 % | ≤ 3 % |
| P6 Sc z² offsets (spreads) | 0.60–0.87 Å (2.5–3.8 Å²) | 0.12–0.31 Å (1.3–2.0) | 0.75–0.97 Å (2.8–3.4) | 0.18–0.23 Å (1.2–1.4) | descriptive |
| P6 other host WFs, max offset | 0.27 Å | 0.23 Å | 1.04 Å | 0.19 Å | descriptive |
| Ω total | 206.6 Å² | 164.0 Å² | 317.6 Å² | 227.5 Å² | |

- **Decision rule: W2 for both hold-outs.** Both variants pass P1–P5 in both ribbons.
- **H4 (W1 vs W2, 601-point Γ–X̄ path, over [VBM − 0.5, CBM + 0.5], which both variants freeze):**
  `r110_8` 3.56 meV (at k 0.44 of Γ–X̄), `r100_6` 0.84 meV. Both ≤ 5 meV, so nothing goes back to the authors.
- **`r110_8` W2 sits closer to the limits than the calibration ribbons did** (A4b-1 W2: ≤ 3.3 meV, masses ≤ 0.9 %).
  - P2 4.71 of 5 meV, at the upper edge-pair band (band 130) at CBM + 0.46 eV, k 0.43 of Γ–X̄, near the top of
    the device window. W1 has 1.27 meV there.
  - CB mass +2.9 % of 3 %. The edge-pair CB is heavy (m* 0.924 m_e DFT, 0.951 W2), fitted on 3 points of the
    31-point path.
  - These pass as registered; nothing is retuned.
- **W2 again removes the Sc z² delocalisation:** z² offsets 0.12–0.31 Å against W1's 0.60–0.97 Å. The `r110_8`
  edge Sc s WFs sit 1.52 Å out into the vacuum (W1 1.55), as in r110_4.
- **Masses (PBE, DFT path):** `r110_8` VB −0.303, CB (edge pair) 0.924 m_e; `r100_6` VB −0.200, CB 0.205 m_e. Gaps
  0.2848 / 0.6172 eV, all direct at Γ.
- The W2 Hamiltonians `twin/w90/<r>/W2/<r>_W2_hr.dat` are the PBE-Wannier hold-outs for D7 (claims gate (b)).

## Results: A5-3 `r110_8` HSE06 gap point (job 143636, 2 nodes, 9 h 00 min, 1,152 core-h; 2026-11-15)
`parse_a5.py hse_points` (→ `hse/a5_hse.json`): the registered `!!` check, then the gap from the IBZ grid,
`parse_a4b.vacuum_level` on the HSE potential, and `parse_a4b.first_nonmodel` (set B with the edge Sc s) on
`projwfc.x` of the HSE SCF.
- **Converged:** 7 EXX outer steps, dexx 5.3e-10 Ry (`!!` −1801.22217584 Ry).
- **HSE gap 0.921 eV, direct at Γ** (NC PBE twin 0.285; r110_4 HSE 1.06).
- **E_vac − VBM 3.576 eV** (PBE 3.115).
- **Band edges:**
  - VB interior N p (edge weight 0.12 = uniform, N p 0.92);
  - **CB still the Sc-s edge pair** (edge weight 0.89; Sc s 0.66, Sc d 0.19, |ψ|² 0.79), so **LR3 holds in HSE**.
    The two edge bands are 2 meV apart at Γ.
- **Non-model check against the device window (registered: reported, not a stop):**
  - **At s_out ≥ 0.25 (and 0.20), the first non-model state is the CB edge pair itself:** s_out 0.273 in HSE,
    against 0.239 in PBE.
    - The same crossing as r110_4 in A4b-2, where H3 fired on the edge pair (amendment 2).
    - Its s_out is mostly spill-out into the vacuum (1 − |ψ|² = 0.21). The edge Sc s is in the model, and the
      [110] edge blocks carry it (W1, amendment 2).
    - It is not the Sc d + N p interior band the registered limitation describes.
  - **The interior vacuum-reaching state** (first at s_out ≥ 0.30) is at CBM + 1.312 eV, outside the device window
    (PBE twin + 1.105).
  - **To word at the claims gate:** by the registered rule, `r110_8` carries the "non-model state in the device
    window" limitation, and that state is the LR3 edge CB.
- **Cost: 1,152 core-h,** at the top of the registered 450–960 plus 10–20 % for two nodes.

## Results: A5-3 `r100_6` HSE06 gap point (job 143637, one node, 20 h 14 min, 1,295 core-h; 2026-11-16)
`parse_a5.py hse_points r100_6` (→ `hse/a5_hse.json`; the `r110_8` entry is unchanged by the rerun):
- **Converged:** 8 EXX outer steps, dexx 2.4e-10 Ry (`!!` −2702.26096049 Ry), 10:15 on 16 Nov, inside the 24 h limit.
- **HSE gap 1.629 eV, direct at Γ** (NC PBE twin 0.617; HSE − PBE + 1.01 eV). [100] trend: n = 3 1.96 → n = 6 1.63 →
  sheet 1.415 eV.
- **E_vac − VBM 3.906 eV** (PBE 3.351).
- **Band edges, both interior** (uniform share 0.167):
  - VB N p 0.95, edge weight 0.19;
  - CB Sc d 0.97, edge weight 0.07;
  - VB−1 at VBM − 0.115 eV leans to the edges (edge weight 0.40, N p 0.91); CB+1 at CBM + 0.47 eV is interior (0.18).
- **Vacuum-like state against the device window (registered: reported, not a stop): outside.**
  - First non-model state at **CBM + 0.574 eV** (band 195 at Γ; s_out 0.60, |ψ|² 0.70), the same at s_out ≥ 0.20,
    0.25 and 0.30. The window ends at CBM + 0.5, so the margin is 74 meV.
  - HSE lowers it relative to the CBM by 0.14 eV (PBE twin + 0.714), as in n = 3 (+ 0.52 → + 0.35) and the sheet
    (+ 0.885 → + 0.74).
  - So in [100] HSE it is inside the window at n = 3 and outside at n = 6. Where it crosses between them is not
    computed. For the claims gate: [100] SK results narrower than n = 6 may carry the registered limitation; the
    wording is decided there.
- **Cost: 1,295 core-h**, above the registered 480–1,030: the 8 outer steps ran ≈ 2.3 h each against the ≈ 1.1–1.3 h
  of the scaling.

## A5 complete (2026-11-16)
- **A5 actual ≈ 3.84k core-h** (registered 1.9–4.9k). Project spent 9,547 (registered total 7.5–10.6k).
- Inputs to the claims gate (D7, judged there, not here):
  - PBE-Wannier hold-outs: W2 for both (`twin/a4b_criteria.json`, `twin/a5_h4.json`); `r110_8` W2 near the P2 and
    mass limits (4.71 meV, + 2.9 %).
  - HSE gap points: `r110_8` 0.921 eV, `r100_6` 1.629 eV, both direct at Γ (`hse/a5_hse.json`).
  - Limitations by the registered rule: `r110_8`, the LR3 edge CB itself scores s_out 0.273; `r100_6`, none (vacuum
    state 74 meV outside the window).
