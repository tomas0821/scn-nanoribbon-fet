# A4b: the edge Wannier recipe and HSE06-W90 on the four calibration ribbons

**Status: FROZEN 2026-11-13 (approved), before any A4b output existed.** Any later change is a
dated, disclosed amendment below the criteria, never an edit of them. A4b-1 was submitted the same day (Run). The recipe was designed **after** seeing the A4a PBE (PAW) projections, so
Methods must say that. The criteria are fixed before any Wannier output exists.

Ribbons (decisions of 12–13 Nov): [100] n = 3 bare (15 Å), [100] n = 3 H (20 Å), [110] n = 4 bare (20 Å),
[110] n = 4 H (20 Å). Geometries are the A4a PAW relaxations.

## What the recipe has to handle (A4a Gate M, PAW PBE; `dft/A4a/edge_bands.py`)
| Ribbon | Band edges | Extra orbitals needed | First non-model empty state |
|---|---|---|---|
| [100] n = 3 bare | VB N p, CB Sc d; no edge state at either | none | **CBM + 0.52 eV** (Γ; \|ψ\|² 0.76) |
| [100] n = 3 H | VB top Sc–H band (H s 67 %), CB Sc d | H s | CBM + 0.90 eV (\|ψ\|² 0.84) |
| [110] n = 4 bare | VB N p, **CB = Sc-s edge pair** (Sc s 53–66 %, \|ψ\|² 0.82) | Sc s on the 2 edge Sc | CBM + 1.22 eV (\|ψ\|² 0.76) |
| [110] n = 4 H | VB top Sc–H band (H s 68 %), CB Sc d | H s | CBM + 0.76 eV (\|ψ\|² 0.82) |

- **"Non-model state"** means an empty state whose weight lies largely outside the projected basis. The measure is
  s_out = \|ψ\|² × (Löwdin share on orbitals outside the projection set) + (1 − \|ψ\|²), and a state counts when
  s_out ≥ 0.25. The PAW data split cleanly: the vacuum-like states score 0.33–0.41 and the [110] edge pair 0.21–0.24
  (with edge Sc s in the set). **Any threshold in 0.25–0.30 gives the same windows**; at 0.20 the [110] edge pair
  would count as non-model.
- These vacuum-like states (Sc s ≈ 0.3, 16–24 % outside the atomic spheres, spread over the ribbon faces, not the
  edges) are the likely cause of the sheet's z² delocalisation (A2: frozen top only 0.5 eV below the vacuum level,
  so the z² bands were disentangled together with vacuum-like states). In the ribbons they start only
  **0.5–0.9 eV above the CBM**.

## Recipe
1. **Code path:**
   - **Pseudopotentials:** the NC set (SG15 Sc, Dojo N, SG15 H; 80/320 Ry), at the PAW-relaxed geometries (the
     A2 convention).
   - **PBE twin (A4b-1):**
     - SCF at k 15;
     - NSCF on the full uniform grid (**12 k for [100], 16 k for [110]**, nosym), `pw2wannier90.x`, `wannier90.x`;
     - bands on the 31-k Γ–X̄ path, with `projwfc.x` on them;
     - `pp.x` bare + Hartree potential for the vacuum level.
   - **HSE06 (A4b-2):** fixed occupations, IBZ SCF + `open_grid.x` (route O, A1), q 6 ([100]) / 8 ([110]),
     `ecutfock` 160.
2. **Projections, in the global frame** (x along the ribbon). The SK side applies `Structure.frame` (R45 for
   [110]), as for the sheet.

   | Ribbon | Projections | num_wann |
   |---|---|---|
   | [100] n = 3 bare | Sc d, N p | 96 |
   | [100] n = 3 H | Sc d, N p, **H s** (all 8 H) | 104 |
   | [110] n = 4 bare | Sc d, N p, **Sc s on the 2 edge Sc** | 66 |
   | [110] n = 4 H | Sc d, N p, **H s** (all 4 H) | 68 |

   H s goes on every H, not only on the Sc-bonded ones that carry the VB-top band. Every atom then keeps its valence
   orbitals, which is what an SK H basis would need.
3. **Excluded bands:** the lowest 5 per f.u. (Sc 3s 3p, N 2s), as A2. In the H ribbons these carry 16–18 % H s from
   the N–H bonds (PAW); that is reported, not fixed.
4. **Outer window:** VBM + 10.23 eV (PBE) / + 11.05 eV (HSE), the T2 ribbon rule unchanged. nbnd is the smallest
   count whose band minimum clears it (the survey: 265 / — / — / 201 at L; the 20 Å counts of [100] H and [110] bare
   come from A4b-1's NSCF grid).
5. **Projection-only, `num_iter = 0`**, as A1/A2. That keeps one WF per atomic orbital, so the orbital-resolved edge
   blocks can be read straight from H(R).
6. **Frozen window: two variants, decided by the criteria below.** All valence bands inside the window are frozen
   (no `dis_froz_min`).
   - **W1 (A2's rule exactly as run, `make_a2_w90.py`):** frozen top = min(vacuum level − 0.5 eV, minimum of the
     first band beyond the model count − 0.1 eV). A2 defined "non-model" by band index only. On the sheet the
     vacuum part bound (the first band beyond the count sat 2 eV higher), and it will bind in the ribbons too.
   - **W2 (W1 plus a character test, new for the ribbons):** frozen top = min(W1, first non-model state − 0.1 eV,
     A2's margin), with "non-model" judged from the projections (definition above). From the PAW numbers that is
     about CBM + 0.42 / 0.80 / 1.12 / 0.66 eV.
   - **W3 (diagnostic only):** W1 on projection set S, which adds Sc s on every Sc. It tests whether the z²
     delocalisation goes away when the vacuum-like states get their own WFs. Its basis is not the SK basis.
   - The vacuum level is the potential averaged over the grid points farthest from all atoms (distance ≥ 90 % of the
     cell's maximum), with its spread reported.

## Pre-registered criteria, PBE twin (per ribbon, per variant)
Reference: the NC PBE bands on the 31-k path. Device window = [VBM − 0.5 eV, CBM + 0.5 eV].

| # | Criterion | Limit | Why this limit |
|---|---|---|---|
| P1 | Windows: frozen top ≥ 0.5 eV below the vacuum level; for W2 also ≥ 0.1 eV below the first non-model state; frozen band count ≤ num_wann ≤ outer band count at every k; outer window not capped below VBM + 10.23 eV | as stated | A2 criterion 3; T2 |
| P2 | Max \|E_W90 − E_DFT\| over every band and path point inside the device window | **≤ 5 meV** | ≈ kT/5; the device only sees this range |
| P3 | RMS inside the frozen window | **≤ 3 meV** | A2 criterion 4 |
| P4 | Band gap (path) | **≤ 3 meV** | |
| P5 | Band-edge masses (parabolic fit to W90 and DFT on the same path points within 30 meV of each extremum, at least 3 points) | **±3 %** | keeps the Wannier error well inside the plan's ±10 % mass test (§3 D7) |
| P6 | WF centres and spreads (descriptive): host non-z² WFs vs a 0.10 Å reference; Sc z² offsets; edge Sc s and H s offsets | reported | A2's centre criterion failed and was reported; this one informs, it does not block |

**Decision rule (per ribbon):**
1. W2 passes P1–P5 → **W2**. It uses the SK-compatible basis and freezes no vacuum-like state.
2. Otherwise, W1 passes P1–P5 → **W1**. The edge blocks then carry the z² vacuum character, and that is disclosed.
3. Otherwise → the plan's LR6 route (recount bands, narrow the window). It goes back to the authors first.
4. W3 is reported only. If it alone passes, using an Sc s device basis is a separate decision.

## HSE stage (A4b-2): criteria registered now
- **H1 Windows:** the same rule on the HSE eigenvalues. The vacuum level comes from the HSE potential (`pp.x`, as
  A2); non-model states come from `projwfc.x` on the HSE SCF. That should work on EXX outputs, which only need the
  saved wavefunctions, but it is untested here. Fallback: match by band index at the same k from the PBE twin.
- **H2:** W90 reproduces the HSE grid eigenvalues inside the frozen window to **≤ 1 meV** (A2 criterion 5).
- **H3 (physics flag; stop if triggered):** HSE pushes the Sc d bands up more than s-like vacuum states, so the
  vacuum-like states (PBE: CBM + 0.5–0.9 eV) may sit inside the HSE device window or below the CBM. **If any HSE
  non-model state lies inside the device window, A4b stops for a decision.** The Sc d + N p (+ edge) basis would
  not carry the device CB there; the options would be an Sc s device basis, or reporting it as a limitation.

## What this does not settle
- **The z² delocalisation, if W2 fails.**
- **How the extra orbitals enter the SK interpolation.** The edge Sc s and H s would be edge-local orbitals in the SK
  ribbon builder. The DFT-ribbon device Hamiltonians use the Wannier H directly, so they are unaffected. The SK H
  basis is still an open item.
- **The 25 Å question.** Settled 13 Nov after the A4b-1 decision: 20 Å stays (decision, option 1); the 8–12 meV
  edge-band residual goes into the error budget (`dft/A4a/README.md` → "25 Å check").

## Cost (`cost_a4.py` law, NC 80 Ry; staging will fix it)
- **A4b-1 (PBE twin + recipe test):** ≈ 170–330 core-h for the four ribbons. Breakdown:
  - [100] bare 35–70;
  - [100] H 67–134;
  - [110] bare 32–64;
  - [110] H 32–65.
- About half of that is the PBE-W90 the budget already carries. **≈ 80–160 is new** (the SCFs, band paths and
  projections that set the windows). This corrects the "≈ 90–190" quoted earlier, which covered only the
  NSCF + `pw2wannier90.x` step.
- Project total ≈ 4.9–12.4k core-h.
- A4b-2 (HSE) as in `cost_a4.py`: **0.9–3.0k** for these four ribbons (HSE-W90 plus PBE-W90: [100] bare 186–629,
  [100] H 351–1,161, [110] bare 153–516, [110] H 233–651).

## Run (A4b-1; inputs generated 2026-11-13, not uploaded)
```bash
~/.venvs/scnfet/bin/python -P dft/A4b/make_a4b.py     # inputs/, w90/<r>/<r>_{B,S}.{win,pw2wan}, wf_map.json, tags_<r>.txt
# copy dft/A4b/ to /path/to/scn_fet/A4b/, then, one job per ribbon, side by side:
sbatch -J a4b1_<r> run_a4b.sh <r>                     # scf, pot, bands, projwfc, nscf, w90 -pp, pw2wannier90 (B, S)
# copy back outputs/ and pot_<r>.cube.gz, then locally:
PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A4b/parse_a4b.py windows    # → a4b_windows.json, w90/<r>/W{1,2,3}/*.win
# upload the variant .win files, then:
sbatch run_a4b_w90.sh r100_3 r100_3h r110_4 r110_4h  # the 12 Wannier90 runs, side by side
# copy back the variant .wout/_hr.dat/_centres.xyz; parse_a4b.py criteria (written while job 1 runs)
```
- **Checks done:**
  - every A4b geometry is identical to its A4a input (`nm_r100_3`, `vac20_r100_3h`, `vac20_r110_4`, `nm_r110_4h`);
  - `run_a4b.sh` and `run_a4b_w90.sh` were dry-run with stub binaries (call order, the gzipped cube, the pdos
    clean-up, the S → B links, a rerun skips everything);
  - `parse_a4b.py windows` was tested on the A4a PAW path and projections of `r110_4h` with a synthetic potential
    cube. It recovers the cube's vacuum level exactly, puts the first non-model state at CBM + 0.758 eV (as the
    PAW analysis), and passes the band-count checks.
- **The 25 Å check runs alongside** (`make_a4a.py vac25` → `vac25_r100_3h`, `vac25_r110_4`; `run_a4a.sh
  tags_vac25.txt`, about 27–54 core-h). `parse_a4a.py stage2` reports 20 → 25 Å with the same 5 meV limit.

## Results: A4b-1 job 1 (jobs 142945–142948, 13 Nov; 157 core-h, estimate 170–330)
Every tag finished ([done]); 30–55 min per ribbon on one node. `parse_a4b.py windows` → `a4b_windows.json`, the 12
variant `.win` files (uploaded, md5 match). NC PBE, 31-k Γ–X̄ path.

| Ribbon (vacuum) | Gap (eV) | E_vac − VBM (spread) | First non-model state (s_out ≥ 0.25) | W1 / W3 frozen top | W2 frozen top | N(L) |
|---|---|---|---|---|---|---|
| [100] n = 3 bare (15 Å) | 0.915 | 3.37 eV (5.5 meV) | CBM + 0.522 (band 99, \|ψ\|² 0.77) | CBM + 1.96 (vacuum) | CBM + 0.42 | 259 |
| [100] n = 3 H (20 Å) | 0.777 | 5.80 eV (14.9 meV) | CBM + 0.923 (105, 0.85) | CBM + 4.52 (vacuum) | CBM + 0.82 | 261 |
| [110] n = 4 bare (20 Å) | 0.419 | 3.26 eV (23.3 meV) | CBM + 1.227 (68, 0.77) | CBM + 2.34 (vacuum) | CBM + 1.13 | 207 |
| [110] n = 4 H (20 Å) | 0.287 | 4.64 eV (55.4 meV) | CBM + 0.761 (68, 0.82) | CBM + 3.85 (vacuum) | CBM + 0.66 | 199 |

- **The NC non-model states sit where the PAW analysis put them** (PAW 0.52 / 0.90 / 1.22 / 0.76 eV above the CBM;
  NC within 0.03 eV), so the W2 tops match the frozen-recipe estimates (0.42 / 0.80 / 1.12 / 0.66).
- **Threshold:** 0.25 and 0.30 give identical windows. At 0.20 only the [110] bare first non-model state moves, onto
  the CB edge pair (as stated in the recipe).
- **W1 is bound by the vacuum level in all four** (the first band beyond the model count lies higher).
- P1 band counts pass for every variant; the outer window (VBM + 10.23 eV) is not capped.
- **H raises E_vac − VBM** by 2.4 eV ([100]) and 1.4 eV ([110]). It only sets the W1/W3 tops here. Edge dipoles
  are a plausible cause; not tested.
- **Path coverage:** for `r100_3h` the band path (150 bands) reaches VBM + 5.22 eV, short of the W1/W3 frozen top
  + 0.1 eV (5.40). Its W1/W3 P3 therefore covers the frozen window only up to VBM + 5.22 eV. That matters only if
  W2 fails for this ribbon.
- **N(L)** (VBM + 10.28 eV, NSCF grid) is the 20 Å band survey the A4b-2 nbnd needs. [110] bare: L + 1 eV lies
  beyond its 225 NSCF bands.

## Results: A4b-1 job 2, the Wannier90 variants (job 143002, 13 Nov; cn004, 17 min, 18 core-h)
All 12 runs ended "All done". `parse_a4b.py criteria` → `a4b_criteria.json`. Limits: P2 ≤ 5 meV, P3 ≤ 3 meV,
P4 ≤ 3 meV, P5 ±3 %.

**Decision (by the frozen rule): W2 for all four ribbons.** W2 passes P1–P5 everywhere, with margin.

| Ribbon | Variant | P2 device max | P3 frozen RMS (max) | P4 gap | P5 m* VB / CB | Verdict |
|---|---|---|---|---|---|---|
| [100] n = 3 bare | **W2** | 3.29 meV | 0.10 (2.0) meV | −0.06 meV | +0.1 / +0.9 % | **pass → used** |
| | W1 | 0.23 | 1.63 (45.3) | −0.06 | +0.1 / +0.3 % | pass |
| [100] n = 3 H | **W2** | 0.27 | 0.17 (4.1) | +0.03 | 0.0 / +0.2 % | **pass → used** |
| | W1 | 1.15 | 1.49 (47.8) | +0.03 | 0.0 / +0.3 % | pass |
| [110] n = 4 bare | **W2** | 2.13 | 0.54 (7.5) | −0.06 | −0.6 / −0.4 % | **pass → used** |
| | W1 | 1.25 | 0.47 (10.1) | −0.05 | −0.1 / +1.4 % | pass |
| [110] n = 4 H | **W2** | 0.23 | 0.23 (5.7) | −0.03 | +0.4 / +0.4 % | **pass → used** |
| | W1 | 0.30 | 1.43 (20.4) | −0.04 | +0.4 / +0.1 % | pass |

DFT masses at the edges (m₀): [100] bare VB −0.21, CB 0.28; [100] H VB −2.15, CB 0.27; [110] bare VB −27.1 (very
flat), CB 0.83; [110] H VB −0.95, CB 0.33. P5 holds even on the flat [110] VB (the watch item did not bite).

**P6, centres and spreads (descriptive):**

| Ribbon | Sc z² offset, W2 (W1) | z² spread, W2 (W1) | Host non-z² max offset, W2 (W1); count > 0.10 Å | Ω total, W2 (W1) |
|---|---|---|---|---|
| [100] bare | 0.18–0.22 Å (0.75–0.91) | 1.2–1.4 Å² (2.8–3.2) | 0.18 Å (1.01); 4 of 84 | 114 (165) Å² |
| [100] H | 0.11–0.19 (0.26–0.60) | 1.2–1.3 (2.5–3.8) | 0.19 (0.46); 12 | 134 (237) |
| [110] bare | 0.13–0.32 (0.60–0.81) | 1.3–2.0 (2.3–3.9) | 0.26 (0.28); 2 | 89 (108) |
| [110] H | 0.12–0.17 (0.35–0.91) | 1.1–1.3 (2.8–5.0) | 0.16 (0.28); 6 | 83 (154) |

- **W2 largely removes the z² delocalisation.** The Sc z² WFs move from 0.3–0.9 Å off their atoms (W1, A2's rule)
  to 0.1–0.3 Å, and their spreads halve. Supporting the recipe's hypothesis: the delocalisation came from freezing
  the vacuum-like states together with the z² bands. It is still above A2's 0.10 Å reference, so the centre test
  stays descriptive, as registered.
- **Edge Sc s ([110] bare):** both WFs sit **1.52 Å outside their Sc atom**, outward in y (−1.22 Å) and out of the
  plane (−0.91 Å), spread 3.4 Å². This is the vacuum-reaching CB edge state (18 % spill-out in Gate M). An SK
  "edge Sc s orbital" would have to carry that off-site centre (open item: the SK H basis).
- **H s:** the two kinds of H differ.
  - Sc-bonded (Sc–H 1.86–1.88 Å): compact, 0.10–0.17 Å off the H, spread 1.4–1.5 Å².
  - N-bonded (N–H 1.02–1.04 Å): spread 2.8–3.4 Å². The centre sits 0.56–0.73 Å **beyond** the H, on the N–H axis
    pointing away from N (perpendicular part ≤ 0.05 Å), i.e. into the vacuum. The N–H bonding state lies partly in
    the excluded bands (16–18 % H s there, recipe item 3); whether that explains the outward shift is not tested.
- **W3 (Sc s on every Sc; diagnostic):** it does not localise the z² WFs better than W2 does. The interior Sc s
  WFs are diffuse (spreads up to 56–81 Å²: vacuum-like states). In both H ribbons disentanglement hit 4000
  iterations without converging. [110] H W3 fails P2 (11.0 meV), P3 (20.1 meV RMS) and P5 (VB +28 %). No reason to
  propose an Sc s device basis.

**What it does not show:**
- whether HSE keeps the vacuum-like states out of the device window (A4b-2, H3);
- W2's frozen top is the character test at 0.25. The 0.20–0.30 robustness was checked on the windows, not on
  the Wannier results;
- the vacuum residual (8–12 meV on the edge bands at 20 Å; `dft/A4a/README.md` → "25 Å check") is a DFT
  convergence error, separate from P2.

**A4b-1 total: 175 core-h** (estimate 170–330).

## A4b-2 staging (HSE06 + Wannier90; inputs generated 2026-11-13; **submitted 21:13 the same day, approved**: jobs 143006–143009)
Rests on: recipe **W2** for all four ribbons (A4b-1 decision) and **20 Å** (decision, 13 Nov). Criteria H1–H3 as
registered above; nothing in them changes.

**Settings** (`make_a4b.py hse` → `inputs/{hse,hpot,hproj,hopen}_<r>.in`, `w90/<r>/H/<r>_H.{win,pw2wan}`,
`tags_hse_<r>.txt`):
- **Method, as A1/A2:** HSE06, NC 80/320 Ry, `ecutfock` 160, fixed occupations, default divergence treatment. SCF on
  the IBZ of the Wannier grid (7 of 12 k for [100], 9 of 16 for [110]), q 6 / 8 along x, then `open_grid.x` to the
  full grid.
- **SCF parameters:** `conv_thr` 1e-9, local-TF mixing β 0.3 (as A4b-1), `diago_full_acc`.
- **Same structure as A4b-1:** geometries identical to the PBE twin (checked), and the projection set is B (checked
  against `<r>_B.win`).
- **nbnd 267 / 273 / 214 / 205** (`parse_a4b.py hsurvey`): the smallest count whose minimum on the A4b-1 PBE NSCF grid
  lies at or above VBM + 10.23 + 0.5 eV.
  - **Why the margin:** the T2 rule assumes the band minima near the window rise by +0.82 eV from PBE to HSE. On the
    A2 sheet they rise by +0.69 to +1.40 eV. A vacuum-like state lifted less could fall inside VBM + 11.05 eV and cap
    the outer window. That is an H1 failure and means a full rerun.
  - **What it costs:** +3–9 % EXX (nbnd²). Without the margin the counts would be 259 / 261 / 207 / 199.
- **Parallel layout:** 63 ranks, one pool per IBZ k (7 × 9 / 9 × 7). With 8 pools, [110]'s ninth k would run
  alone and double the wall time.
- **Memory:** the EXX buffers are ≈ 17–30 GB per pool copy; the nodes have 1 TB.

**Jobs:**
- **Job 1, one per ribbon** (`run_a4b2.sh <r>`, 24 h on one node). Tags, in order:
  - `hse` (pw.x);
  - `hpot` (pp.x, the HSE potential → E_vac);
  - `hproj` (projwfc.x on the HSE SCF → non-model states, H1/H3);
  - `hopen` (open_grid.x);
  - `hw90pp` (fills the open_grid k list into the `.win`, refuses a wrong count, then `-pp`);
  - `hp2w` (pw2wannier90).
  - Restartable. An `hse` rerun wipes everything downstream. There is no EXX restart: `max_seconds` 80,000 s
    (22.2 h) is only a guard.
- **Laptop:** `parse_a4b.py hwindows` → H1, the H3 stop flag, the W1/W2 `.win` files.
- **Job 2** (`run_a4b2_w90.sh`, 2 h): W2 (the recipe) and W1 (reported).
- **Laptop:** `parse_a4b.py hcriteria`. It checks H2 and reports the HSE gap, edge masses and centres
  (descriptive).

**Cost** (`cost_a4.py`, staged nbnd): HSE06 SCFs **1,028–2,528 core-h**, plus ≈ 30–80 for the other steps.

| Ribbon | nbnd | HSE SCF (core-h) | Wall on one node |
|---|---|---|---|
| [100] n = 3 bare | 267 | 206–506 | 3.2–7.9 h |
| [100] n = 3 H | 273 | 356–876 | 5.6–13.7 h |
| [110] n = 4 bare | 214 | 243–598 | 3.8–9.3 h |
| [110] n = 4 H | 205 | 223–548 | 3.5–8.6 h |

The earlier quote, 0.9–3.0k, included the PBE-W90 that A4b-1 has now done. With the staged band counts, the
project total is **5.2–12.0k core-h** (spent 1,667).

**Checks done:**
- **Inputs:** the 45 A4b-1 input files regenerate byte-identical. Re-parsing A4b-1 gives identical windows and
  criteria JSON and identical variant `.win` files.
- **Dry runs with stub binaries:**
  - call order, ranks and pools (63 / `-nk 7` or `-nk 9`);
  - the cube is gzipped and the pdos files are removed;
  - the k list is filled (12 / 16 points); an 11-point list is refused and the placeholder kept;
  - an `hse` rerun redoes every downstream step; a finished job skips everything;
  - the Wannier90 job links the HSE `.amn/.mmn/.eig` into W1/W2.
- **Proxy test of `hwindows` / `hcriteria` on the PBE twin data:**
  - the frozen windows equal A4b-1's, and the outer window moves up by the intended 0.82 eV;
  - H2 gives 0.06–0.26 meV;
  - the Wannier edge masses equal the DFT ones;
  - the cap flag fires when the proxy has too few bands (the [110] bare NSCF's 225).

**Reading the results:**
- At s_out 0.20 the [110] bare CB edge pair itself counts as non-model, so H3 "at 0.20" will read "stop" for that
  ribbon. The registered threshold is 0.25; 0.20 and 0.30 are reported only.
- **Untested:**
  - `projwfc.x` on an EXX run. If it fails, H1's registered fallback applies (band index from the PBE twin), and it
    goes back for a decision;
  - `open_grid.x` on a 1D grid (A2 used it in 2D).

**Run (after approval):**
```bash
# upload inputs/{hse,hpot,hproj,hopen}_*.in, tags_hse_*.txt, run_a4b2*.sh, w90/<r>/H/<r>_H.{win,pw2wan}; then
sbatch -J a4b2_<r> run_a4b2.sh <r>                    # four jobs, side by side
# copy back outputs/{hse,hpot,hproj,hopen,hw90pp,hp2w}_*.out, hpot_<r>.cube.gz, w90/<r>/H/<r>_H.{win,eig}; then
PYTHONPATH=. ~/.venvs/scnfet/bin/python -P dft/A4b/parse_a4b.py hwindows      # H1, H3 → a4b2_windows.json
# upload w90/<r>/H/W?/*.win; sbatch run_a4b2_w90.sh r100_3 r100_3h r110_4 r110_4h; copy back; then hcriteria
```

## Results: A4b-2, first two ribbons (2026-11-14) — **H3 stop flag on [100] bare**
- **Done:** 143006 r100_3 (7.1 h, 454 core-h; 5 EXX cycles, dexx 7.7e-10 Ry) and 143009 r110_4h (13.5 h, 861 core-h;
  5 cycles, 7.4e-10).
- **Running at 13:06:** 143008 r110_4 (sixth cycle) and 143007 r100_3h (about two cycles left; near the 80k s
  `max_seconds` guard). `hwindows` now refuses an output without QE's final "!!" energy, because a max_seconds stop
  may still print JOB DONE.

| Ribbon | HSE gap (IBZ grid) | PBE twin gap | E_vac − VBM | First non-model (s_out ≥ 0.25) | H3 | W2 frozen top | H1 |
|---|---|---|---|---|---|---|---|
| [100] n = 3 bare | 1.964 eV | 0.915 eV | 3.93 eV | **CBM + 0.351 eV** (Γ, band 98, \|ψ\|² 0.70) | **STOP** | CBM + 0.251 | ok |
| [110] n = 4 H | 1.487 eV | 0.287 eV | 5.41 eV | CBM + 0.642 (band 68, \|ψ\|² 0.78) | clear | CBM + 0.542 | ok |

- **[100] bare, the state that stops it:** Sc d 54 %, Sc s 35 %, N s 11 %, 30 % outside the atomic spheres. It lies
  on the faces (edge share 0.15), not the edges. It disperses (+0.395 eV at k = 1/12), and it is the PBE twin's band
  99 at Γ (CBM + 0.522).
  - HSE raised the Sc d CB by 1.05 eV but this state by only 0.88 eV: the case H3 anticipated.
  - s_out 0.62 is far above any threshold in 0.20–0.30.
- **Outer window:** not capped in either ribbon (top band minima VBM + 11.59 / 11.68 eV); the 0.5 eV nbnd margin held.
- **Cost:** the law fits [100] bare (454 in 206–506) but [110] H ran ×1.6 over the top of its range.
- Per the registered rule, A4b stops for a decision; no Wannier90 job submitted.

## Amendment 2026-11-14 (decision after the H3 stop; written before any W3 HSE output exists)
The frozen criteria above are unchanged. This adds what happens after H3 fires.

1. **Device basis when H3 fires.** The device Hamiltonian of that ribbon uses **W3**: projection set S (set B plus
   Sc s on every Sc), with the W1 window rule (frozen top = min(E_vac − 0.5 eV, first band beyond the set-S count −
   0.1 eV)). It must pass:
   - H1 (counts, vacuum margin, outer window not capped);
   - H2 (≤ 1 meV on the HSE grid inside its frozen window);
   - and its frozen window must contain the whole device window (VBM − 0.5 to CBM + 0.5 eV).

   If W3 fails any of these, the ribbon falls back to option 1: W2, plus the band reported as a limitation with a
   top-of-barrier estimate of its share. There is no further tuning.
2. **W2 is still run for that ribbon** and supplies the orbital-resolved blocks for the SK side (the SK basis has no
   Sc s). That the SK interpolation does not carry this band is stated as a limitation.
3. **Applies now to [100] n = 3 bare.**
4. **Pre-registered for the A5 hold-outs** ([110] n = 8, [100] n = 6), before any of their output exists: the same
   H1–H3 checks on their HSE runs, and the same rule (items 1–2) if H3 fires.
5. **Where the band sits does not change.** The state is a property of the free-standing ribbon. A gate dielectric
   would likely push such vacuum-reaching states up; the manuscript says so and treats the free-standing result as
   the conservative case.
6. **Sheet check (diagnostic):** `projwfc.x` on the A2 HSE and PBE runs, to see whether the M-ScN sheet has the same
   kind of state near its CBM. This is reported only; it changes no A2 or gate (a) result.

## Results: A4b-2 Wannier90 for [100] bare and [110] H (job 143231, 14 Nov; 23 min) and the sheet check (job 143232)
`parse_a4b.py hcriteria` → `a4b2_criteria.json`; `dft/A2/sheet_nonmodel.py` → `dft/A2/a2_sheet_nonmodel.json`.

| Ribbon | Variant | H1 | H2 (HSE grid, frozen window) | Gap (W90 path) | m* VB / CB | Sc z² offset | Ω | Use |
|---|---|---|---|---|---|---|---|---|
| [100] n = 3 bare | **W3** (set S, frozen to CBM + 1.47 eV) | ok, device window inside | **0.275 meV** (290 values) | 1.964 eV, direct Γ | −0.194 / 0.273 | 0.20–0.34 Å | 208 Å² | **device Hamiltonian** (amendment) |
| | W2 (frozen to CBM + 0.25) | ok | 0.279 meV | 1.964 | −0.194 / 0.276 | 0.31–0.41 Å | 133 | **SK blocks** |
| | W1 | ok | 0.270 meV | 1.964 | −0.194 / 0.274 | 0.73–0.92 Å | 173 | reported |
| [110] n = 4 H | **W2** (frozen to CBM + 0.54) | ok | **0.058 meV** (236 values) | 1.487 eV, indirect (VBM X̄, CBM Γ) | −0.863 / 0.317 | 0.22–0.27 Å | 95 | **device + SK** |
| | W1 | ok | 0.054 meV | 1.487 | −0.862 / 0.318 | 0.39–1.05 Å | 183 | reported |

- **Convergence:** [100] bare W3 ended at the 4000-iteration cap with a last change of 2e-10 Å² (tolerance 1e-10),
  so it is converged in practice. [110] H W1 hit the cap at 7e-7 (W1 is reported only).
- **What W3 adds:** it reproduces the HSE eigenvalues to 0.28 meV across a frozen window that contains the
  nearly-free-electron band at CBM + 0.35 eV. The device Hamiltonian of [100] bare now carries that band.
- **The edge masses are close to the PBE twin's DFT masses:** [100] bare VB −0.209, CB 0.277; [110] H −0.95, 0.33.
- **The sheet has the same kind of state** (diagnostic, amendment item 6):
  - HSE: CBM + 0.743 eV at Γ, band 18, \|ψ\|² 0.71, Sc d 57 %, Sc s 33 %, N s 10.5 %;
  - PBE: CBM + 0.885 eV;
  - so it lies outside the sheet's device window.
  - **But it lies inside A2's frozen window (CBM + 1.8 eV) and inside the ±1 eV window that gate (a) scored.** The
    sheet Wannier bands the SK fit targeted therefore contain this band, which an Sc d + N p SK model is not built
    to carry. That was not known when gate (a) ran; whether it affected the fit is not checked. Gate (a) stands as
    run.
- **Trend:** the state sits at CBM + 0.35 eV in the narrow bare [100] ribbon and + 0.74 eV in the sheet; [110] H
  + 0.64. So the wider [100] n = 6 hold-out probably lies between, and H3 may or may not fire there; the amendment
  covers it either way.

## Results: A4b-2 [110] n = 4 bare (job 143008, 16 h 50 min, 1,077 core-h; Wannier90 job 143234, 19 min)
- **HSE:** gap 1.058 eV on the IBZ grid (PBE twin 0.419); 6 EXX cycles, final dexx 9.6e-11 Ry. Outer window not capped
  (top band VBM + 11.33 eV).
- **H3 fired, on the CB edge pair itself.** The pair is a modelled state (edge Sc s 62–66 %; only 5–9 % of its
  projected weight lies outside set B). Its s_out rose from 0.23–0.25 (PBE) to 0.26–0.28 (HSE) because \|ψ\|² fell
  from 0.82 to 0.79 (more spill into the vacuum). At s_out 0.30 the first non-model state is CBM + 1.446 eV, so no
  vacuum-like band lies in the device window. By the frozen rule H3 fires all the same, and the amendment's W3 was run.

| Variant | H1 | H2 (grid) | Gap, 601-point W90 path | CB min | m* VB / CB | z² offset | Ω | Disentanglement |
|---|---|---|---|---|---|---|---|---|
| W1 (set B, frozen CBM + 2.16) | ok | 0.051 meV | 1.053 eV | Γ | −19.9 / 0.95 | 0.60–0.76 Å | 117 Å² | converged (1467 it.) |
| W2 (frozen CBM − 0.10) | ok | 0.052 meV | 1.503 eV (misses the edge CB) | Γ | −20.0 / 1.08 | 0.26–0.40 | 101 | converged |
| W3 (set S, frozen CBM + 2.16) | ok, device window inside | 0.052 meV | **0.996 eV** | **k = 0.07 of Γ–X̄** | −19.2 / **0.18** | 0.05–0.44 | **868** | **cap (4000), change 3.5e-9** |

- **Cross-variant check on a 601-point path, device window:**
  - [110] bare **W1 vs W3: up to 222 meV**, and W3's CB dips 57 meV below the DFT CBM between grid points;
  - [100] bare W1 vs W3: 1.5 meV; [100] bare W1 vs W2: 0.7 meV; [110] H W1 vs W2: 1.0 meV.
- **Why W3 fails here:** four interior Sc s WFs (Sc 2, 4, 11, 13) spread to 168–195 Å², fully delocalised; in
  [100] bare W3 the largest Sc s spread is 5.6 Å². H2 tests only the 16 grid points, so it cannot see this.
- **So the amendment's rule picks W3, but W3's interpolation between grid points is wrong.** The PBE twin's W3 for
  this ribbon passed every path criterion, so this is new in HSE. **Back to the authors.** Proposed:
  W1 for the device Hamiltonian and the SK blocks of [110] bare, as a dated deviation, plus a cross-variant check
  pre-registered for the hold-outs.

## Amendment 2026-11-14 (2) (decision after the [110] bare W3 result)
1. **[110] n = 4 bare: W1 is the device Hamiltonian and the source of the SK blocks.** This is a dated deviation
   from amendment (1), which selects W3.
   - Why: W3 passes H1/H2 on the 16 grid points but is wrong between them. It differs from W1 by up to 222 meV in
     the device window, four interior Sc s WFs are delocalised (168–195 Å²), and its disentanglement hit the cap.
   - W1 uses the same frozen window (CBM + 2.16 eV) and set B, which already holds the edge Sc s. It is converged,
     and it agrees with W2 within 0.2 meV over the valence part of the device window (all that W2 freezes).
   - W1's Sc z² WFs sit 0.60–0.76 Å off their atoms (disclosed, as the PBE decision rule's W1 clause says).
   - W2 cannot serve: its frozen top lies below the edge CB.
2. **New check H4, pre-registered for the A5 hold-outs (before any of their output):** every pair of Wannier variants
   run for a ribbon must agree within **5 meV** (P2's limit) over the part of the device window that both freeze, on
   a 601-point Γ–X̄ path. If any pair fails, the ribbon goes back to the authors; no variant is picked automatically.
3. **For the four calibration ribbons H4 is computed and reported** (descriptive; it was not registered before
   their output; `a4b2_criteria.json` → `H4`):
   - [100] bare: W1–W2 0.7, W1–W3 1.5, W2–W3 0.8 meV;
   - [110] H: W1–W2 1.0 meV;
   - [110] bare: W1–W2 0.2, W2–W3 4.6 meV (valence side only), **W1–W3 222 meV** (the case above).

## Results: A4b-2 [100] n = 3 H (job 143007, 21 h 16 min, 1,361 core-h; Wannier90 job 143311, 16 min) and the summary
- **HSE:** gap 1.972 eV (PBE twin 0.777); 5 EXX cycles, final dexx 8.3e-10 Ry. pw.x wall 21.0 h, inside the 22.2 h
  `max_seconds` guard. Outer window not capped.
- **H3 clear:** first non-model state at CBM + 0.829 eV (PBE + 0.923).
- **W2:** H2 0.322 meV, converged (3033 iterations). W1 hit the iteration cap (reported only).
- **H4:** W1–W2 1.5 meV.

**A4b-2 summary: the four HSE06 device Hamiltonians** (`figures/A4b2_HSE_device_Hamiltonians.png`)

| Ribbon | HSE gap (PBE twin) | Device Hamiltonian | H2 | H4 (pairs) | m* VB / CB (W90 path) | Note |
|---|---|---|---|---|---|---|
| [100] n = 3 bare | 1.964 eV (0.915), direct Γ | **W3** (amendment 1) | 0.275 meV | ≤ 1.5 meV | −0.194 / 0.273 | NFE band at CBM + 0.35 eV carried; W2 gives the SK blocks |
| [100] n = 3 H | 1.972 eV (0.777), direct Γ | **W2** | 0.322 meV | 1.5 meV | −2.21 / 0.273 | |
| [110] n = 4 bare | 1.058 eV (0.419) | **W1** (amendment 2) | 0.051 meV | W1–W3 222 meV | −19.9 / 0.95 | H3 fired on the edge pair; W1 also gives the SK blocks |
| [110] n = 4 H | 1.487 eV (0.287), indirect | **W2** | 0.058 meV | 1.0 meV | −0.86 / 0.317 | |

- **Cost:** HSE job-1s 3,754 core-h (estimate 1,028–2,528), plus the Wannier90 jobs and sheet check ≈ 63, so
  A4b-2 ≈ 3.8k.
  - The EXX constant per run: 6.0e-9 at 15 Å, 1.05–1.06e-8 at 20 Å (three nodes).
  - `cost_a4.py` is recalibrated (`C_EXX_RIBBON`, "After A4b-2").
