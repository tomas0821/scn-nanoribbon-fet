# A4a — calibration ribbons and Gate M (PBE)

Plan §4 A4a, §3 D4, §10 LR3/LR5. Written 2026-11-11, **before any A4a output exists**. A4a-1 submitted 2026-11-11
(approved) as jobs 142068 (`relax_r100_3`), 142069 (`relax_r110_4`), 142070 (`relax_r100_3h`).

## Ribbons (`make_a4a.py`, from `scnfet.structures.ribbon`)
| Tag | Ribbon | Atoms | W (edge to edge) | Cell (Å) |
|---|---|---|---|---|
| r100_3 | [100] n = 3, bare | 24 (12 Sc, 12 N) | 1.10 nm | 4.397 × 25.99 × 17.17 |
| r110_4 | [110] n = 4, bare | 16 (8 Sc, 8 N) | 1.09 nm | 3.109 × 25.88 × 17.17 |
| r100_3h | [100] n = 3 + 8 H | 32 | 1.47 nm | 4.397 × 25.99 × 17.17 |

- **The H-terminated ribbon is [100]:** each edge atom has one cut bond along the edge normal, so one H per edge atom
  (N–H 1.02, Sc–H 1.83 Å to start) is unambiguous. The [110] H termination stays an open chemistry decision.
- **Geometry convention (as A2):** the A0 PBE sheet scaled by bulk a_exp / a_PBE: a_p 3.1093 Å, planes 2.1685 Å.
  Internal coordinates are relaxed at that fixed period, and the relaxed structures are reused unchanged by A4b
  (HSE-W90) and the device Hamiltonians.
- **Vacuum:** 15 Å in y and z, checked by `vac20`.

## Settings
- **Functional and pseudopotentials:** PBE, SSSP precision as A0 (Sc PAW psl 0.2.3 spn, N Dojo NC, H SG15 = the SSSP
  H, 80/320), ecutwfc 90 / ecutrho 720 Ry (the A0 rule).
- **Smearing:** gaussian 0.005 Ry for every run, NM and magnetic alike, so their energies compare.
- **k-points:** 15 along x (Γ-centred, 8 irreducible); k-test at 30.
- **Relaxations:** BFGS, forces ≤ 1e-3 Ry/bohr, ΔE ≤ 1e-5 Ry; local-TF mixing, β 0.3.
- **Bare ribbons relax from a ferromagnetic start** (nspin 2; 0.1 × 11 e on Sc, 0.2 × 5 e on N), so an edge moment can
  survive or collapse. The H ribbon relaxes with nspin 1 and gets an FM check afterwards.

## Stages (`run_a4a.sh <tag | tag-list file>`, one 64-rank node per job, restartable)
- **A4a-1:** one job per relaxation, side by side: `sbatch -J a4a_<tag> run_a4a.sh <tag>` for `relax_r100_3`,
  `relax_r110_4`, `relax_r100_3h` (listed in `tags_stage1.txt`).
  - A relax counts as done only when its output has "End final coordinates": one stopped by `max_seconds` also
    prints "JOB DONE".
  - Resubmitting an unfinished relax restarts it (`restart_mode='restart'` from its `.save` dir; the old output is
    kept as `.out.<n>`).
- **A4a-2** (`python3 make_a4a.py stage2` after A4a-1 is copied back; 23 tags; one job per ribbon:
  `sbatch -J a4a2_<ribbon> run_a4a.sh tags_stage2_<ribbon>.txt`). At each relaxed geometry:
  - NM, FM, AFM across the edges (left half ↑, right half ↓) and AFM between the planes (lower ↑, upper ↓) at k 15;
  - NM, FM and AFM-edges at k 30;
  - bands Γ–X̄ (31 k; 28 bands per f.u. + 1 per H, which also surveys A4b's band count) and k-resolved `projwfc`
    for NM;
  - the H ribbon: NM + FM;
  - vacuum check: r110_4 NM at 20 Å.
- The generator was tested on synthetic relaxed outputs, and the job script was dry-run with a stub `srun`: all tags in
  order, a rerun skips them all, and the tmp dirs are kept only where they are reused.

## Cost (`cost_a4.py`; the model is in that file's docstring)
| Run | core-h |
|---|---|
| A4a-1: three relaxations | estimate 330–660; **actual 580** (190 + 110 + 280; jobs 142068–142070, 12 Nov) |
| A4a-2: single points, k-test, bands (with the band survey), vacuum | estimate 390–780; **actual 357** (200 + 97 + 61; jobs 142490–142492, 12 Nov) |
| **A4a total** | **937** (plan 600–1,200) |

## Pre-registered Gate M criteria
Edge atoms per cell: r100_3 has 8 (2 edges × 2 planes × Sc + N); r110_4 has 4.
1. **Magnetism.** ΔE_mag = E_NM − min(E_FM, E_AFM-edges, E_AFM-planes), per edge atom.
   - **Magnetic if ΔE_mag > 5 meV per edge atom at k 30** (LR5 threshold).
   - Reported with the site moments.
   - Magnetic → LR5 (spin-resolved edge term, or H edges only; ≈ 8k budget): for a decision.
   - ≤ 5 meV → the ribbon is treated as non-magnetic for A4b/A5, and ΔE_mag is reported.
2. **k-test.** ΔE_mag changes by ≤ 1 meV per edge atom from k 15 to k 30, and the NM gap by ≤ 10 meV. Otherwise the
   later PBE runs use k 30.
3. **PBE gap.** From the NM bands (31 k on Γ–X̄) and the SCF grid: value, direct or indirect, k positions. It is set
   against the sheet's PBE gap at the same geometry convention (0.404 eV, A2).
   - **Metallic** if a band crosses E_F on the path.
4. **Edge states.** For each band-edge state (VB top, CB bottom), the k-resolved weight on the outermost row of each
   edge (both planes).
   - **Edge > 50 % → the gap is set by an edge state → LR3:** characterise it; the H ribbon becomes the device case
     and both chemistries are reported.
5. **Vacuum.** r110_4 NM gap 15 → 20 Å: ≤ 5 meV (as A2 criterion 2). Otherwise A4b/A5 use 20 Å (≈ +30 % HSE cost).
6. **Edge geometry (descriptive, no threshold).** Maximum displacement from the ideal cut, edge bond lengths, and any
   reconstruction.

Gate M is judged per ribbon from these numbers alone. No ribbon HSE (A4b) starts before it is written down.

## Results (2026-11-12; jobs 142490–142492, 357 core-h; `parse_a4a.py stage2` → `a4a_stage2.json`)
The conventions the criteria leave open are in the `stage 2` header of `parse_a4a.py`. They were fixed before
criteria 2–5 were computed; the criterion-1 energies had already been seen in a job status check.

| Criterion | [100] n = 3 | [110] n = 4 | [100] n = 3 + H |
|---|---|---|---|
| 1. ΔE_mag (meV per edge atom) | −0.11 at k 30; FM, AFM-edges, AFM-planes all → 0 μB | −0.14 at k 30; all three → 0 μB | −0.11 at k 15 (FM → 0 μB) |
| → magnetic (> 5) | **no** | **no** | **no** |
| 2. k-test (ΔΔE_mag; NM gap) | 0.00; 0.1 meV: **pass** | 0.00; 0.2 meV: **pass** | not run (as registered) |
| 3. PBE gap (31 k, Γ–X̄) | **0.912 eV**, direct at Γ | **0.409 eV**, direct at Γ | **0.767 eV**, direct at Γ |
| 4. Edge weight VBM / CBM (uniform share) | 0.45 / 0.32 (0.33) | 0.35 / **0.89** (0.25) | host edge row 0.29 / 0.23 (0.25); **with the H atoms 0.97** / 0.23 |
| → edge state sets the gap (> 0.5) | **no** | **yes: the CB bottom** | host row only: no. **With H: yes, the VB top** (open: interpretation to be decided) |
| 5. Vacuum 15 → 20 Å (r110_4 NM gap) | | **12.3 meV: fail** (limit 5) | |

All runs converged; no metallic ribbon. The sheet's PBE gap is 0.404 eV.

**Characterisation (LR3: "characterise").**
- **[110] n = 4: the CB bottom is an edge-state pair** at Γ, +0.409 and +0.431 eV above the VBM. 82 % of the
  weight sits on the two outermost Sc atoms (41 % each); Sc s 53 %, Sc d 16 %; |ψ|² 0.82, i.e. 18 % of the state lies
  outside the atomic orbitals (spill-out into the vacuum). At X̄ the pair is at +1.71 eV. The first interior-like CB is
  at +1.51 eV (27 % edge), so the pair sits ≈ 1.1 eV inside the interior gap. The VB top is N p, 35 % edge: interior.
- **The vacuum failure is that pair.** Between 15 and 20 Å, VBM-aligned at Γ, the pair moves −12.3 meV, the interior
  CB −1.4 meV, the VB ≤ 0.5 meV.
- **H ribbon: the VB top is a pair of Sc–H states** at Γ (0 and −0.012 eV). 67 % of the weight sits on the four
  Sc-bonded H atoms (H s 67 %, Sc d 12 %, N p 10 %, Sc s 8 %). At X̄ the pair joins a four-fold group at −0.21 eV.
  Bands 96–97 (−0.44, −0.32 eV) are also H-derived. The interior VB top (band 98, N p, 28 % edge) is at −0.267 eV, so
  the interior gap is ≈ 1.03 eV. The CB bottom is interior Sc d on the middle atoms.
- **[100] n = 3 bare: no edge state at either band edge.** VB N p (45 % edge), CB Sc d (32 %). The VB−1 state at
  −0.31 eV is 67 % edge but does not set the gap. **The bare [100] gap (0.912 eV) is larger than the H ribbon's
  (0.767 eV)**: the H termination adds the Sc–H states above the interior VB.

**Band survey for A4b (descriptive, at 15 Å).** Smallest N whose band minimum on the path is ≥ VBM + 10.28 eV (the
sheet's PBE level for the accepted HSE window VBM + 11.05 eV):

| Ribbon | N at VBM + 10.28 eV | at + 10.78 | at + 11.28 | Cost model |
|---|---|---|---|---|
| [100] n = 3 | 265 (22.1 per f.u.) | 272 | 281 | 240–288 |
| [110] n = 4 | 174 (21.75 per f.u.) | 184 | 193 | 160–192 |
| [100] n = 3 + H | 234 (19.5 per Sc, H included) | 243 | 257 | 256–304 |

At 20 Å more vacuum states fall below the same level; the survey does not cover that.

**Gate M, per ribbon:**
- **[100] n = 3 bare:** non-magnetic; k 15 adequate; PBE gap 0.912 eV, direct at Γ; no edge state at the band edges.
- **[110] n = 4 bare:** non-magnetic; k 15 adequate; **PBE gap 0.409 eV, set by an in-gap Sc-s edge state → LR3.**
- **[100] n = 3 + H:** non-magnetic; PBE gap 0.767 eV, direct at Γ; **the VB top is a Sc–H edge state if the H
  atoms count as the edge row** (open).
- **Vacuum (criterion 5) fails.** By the rule above, A4b/A5 use 20 Å. **The "≈ +30 % HSE cost" written above was an
  underestimate:** the larger FFT grid alone makes ribbon HSE ×1.65. With option A, the project total becomes 5.4–14.5k
  core-h instead of 3.7–9.3k, before the extra vacuum bands (`cost_a4.py`, "Vacuum").

**Decisions after Gate M (2026-11-12; made after seeing the results, to be disclosed as such):**
1. **Vacuum:** the [110] ribbons use 20 Å in A4b/A5. The [100] bare and H ribbons get the same 20 Å check (NM gap,
   ≤ 5 meV; `make_a4a.py vac20`, job 142557) and use 20 Å only if they fail.
2. **H counts as the edge row** in criterion 4. The H ribbon's VB top is therefore an edge state (0.97).
3. **LR3:** add a [110] n = 4 H ribbon (`r110_4h`: one H per edge atom on the edge normal, 20 Å vacuum;
   `make_a4a.py h110` then `h110 2`) through relax + Gate M before any [110] HSE. [100] A4b can go ahead.
   The termination waits for approval.

## Vacuum checks on the [100] ribbons (2026-11-13; job 142557, cn016, 37 min, 40 core-h; estimate 25–50)
NM SCF at k 15, the relaxed ribbon centred in the 20 Å cell (`make_a4a.py vac20`), compared with `nm_<ribbon>` at
15 Å (`parse_a4a.py stage2` → `a4a_stage2.json` `c5_vacuum`). Limit 5 meV, as registered.

| Ribbon | Gap 15 → 20 Å (eV) | Δ gap | ΔE_tot | Vacuum for A4b/A5 |
|---|---|---|---|---|
| [100] n = 3 bare | 0.9116 → 0.9116 | **0.0 meV**, pass | −0.0 meV | **15 Å** (also the [100] n = 6 hold-out: same edges) |
| [100] n = 3 H | 0.7672 → 0.7849 | **17.7 meV**, fail | −19.1 meV | **20 Å** |
| [110] n = 4 bare (12 Nov) | 0.4084 → 0.3961 | 12.3 meV, fail | +6.6 meV | 20 Å |

**Which bands move** (Γ, eigenvalues aligned on the median shift of all occupied bands;
`vac_bands.py`):
- **[100] H:** the Sc–H pair at the VB top moves **−17.5 / −17.3 meV**; two lower bands (−0.44, −0.32 eV) −16 / −13 meV.
  The interior VB top (−0.267 eV) moves **+0.2 meV** and the CB (interior Sc d) **+0.4 meV**. So the gap change is the
  edge pair alone, as in [110] (pair −12.0 / −11.9, interior CB −1.1, VB ≤ 0.8 meV). The occupied-band shifts spread
  by 4.8 meV (std), against 0.2 / 0.7 meV in the bare ribbons.
- **[100] bare:** every band within 0.3 meV, except CB+1 (+1.43 eV) at +4.9 meV. The 67 %-edge state at −0.31 eV
  moves −0.3 meV.
- **Pattern:** the vacuum sensitivity belongs to the two edge states that reach into the vacuum: the [110] Sc-s pair
  (18 % outside the atoms) and the Sc–H pair (H s 67 %). Interior bands are converged at 15 Å in all three ribbons.
  Mechanism not tested.

**Is 20 Å enough? Not shown.** The registered test compares 15 with 20 Å only. If the edge shift decays as 1/L^p
(L = vacuum), the 20 → 25 Å change would be 6–11 meV for the H pair and 4–7 meV for the [110] pair (p = 3 … 1). So
20 Å may fail the same 5 meV test one step on. A 25 Å NM SCF on these two ribbons would settle it (estimate
27–54 core-h, `cost_a4.py`); it is not staged.

**Band counts at 20 Å are not surveyed.** The 20 Å SCFs carry 120 / 77 bands, too few for the window. The A4b nbnd for
the H and [110] ribbons needs a PBE band survey at 20 Å first; more vacuum puts more free-electron bands below
VBM + 10.28 eV than the 15 Å counts (234 / 174).

**Cost now** (`cost_a4.py`; spent 1,237 core-h):
- remaining option A with this vacuum set: 3,153–10,106; + [110] n = 4 H at 20 Å (A4a 92–223, HSE-W90 + PBE-W90
  168–559): 3,413–10,888;
- **project total 4,651–12,125** (15 Å everywhere would have been 3,766–9,320);
- if 25 Å turns out to be needed for the edge-state ribbons: project 5,348–14,294.
- The A4a estimate for `r110_4h` is bracketed by the two relaxes it resembles: the law underestimated the [110] bare
  relax ×1.5 and the [100] H relax ×5.1.

**Decisions (2026-11-13):** the 25 Å check waits. The [110] n = 4 H ribbon goes ahead: `relax_r110_4h` is job 142699,
then `make_a4a.py h110 2`. The nominal vacuum is host-to-host: counting the H atoms, the [100] H ribbon's real gap to
its images is 12.8 / 11.8 Å (y / z) at 15 Å and 17.8 / 16.8 Å at 20 Å.

## [110] n = 4 H ribbon (`r110_4h`, LR3; 2026-11-13)
- **Relax** (job 142699, 150 core-h; bracket 53–184): converged in 31 BFGS steps. Sc–H 1.860 Å, N–H 1.024 Å. The H
  atoms tilt out of the plane (Sc-bonded 0.94 Å, N-bonded 0.59 Å); H–H on an edge 2.31 → 3.67 Å. Host edge atoms move
  ≤ 0.17 Å; the edges stay mirror images. W incl. H 1.39 nm; the real vacuum counting H is 17.0 Å (y) / 18.2 Å (z).
- **Step 2** (job 142745): NM, FM start, bands (31 k, **nbnd 316**: the empty bands scaled by the cell cross-section
  for the 20 Å vacuum), projwfc. Gate M criteria 1, 3, 4 as registered for the H ribbon (no k 30, no AFM runs);
  criterion 5 does not apply (already at 20 Å).

### Results (job 142745, cn014, 54 min, 57 core-h; `parse_a4a.py stage2`, `edge_bands.py r110_4h`)
| Criterion | [110] n = 4 + H (20 Å) | [110] n = 4 bare, for comparison |
|---|---|---|
| 1. Magnetic | **No**: FM start → 0 μB; ΔE_mag −0.11 meV per edge atom | No (−0.14) |
| 3. PBE gap (31 k Γ–X̄) | **0.302 eV, indirect** (VBM at X̄, CBM at Γ); direct Γ 0.998, X̄ 0.958; k 15 grid 0.320 (misses X̄) | 0.409 eV, direct Γ |
| 4. Band edges | **CB bottom interior Sc d** (88 %; edge + H 0.14, uniform 0.20; |ψ|² 0.99). **VB top = Sc–H edge band** (H s 68 %, Sc d 24 %; edge + H 0.97) | CB bottom = Sc-s edge pair (0.89; |ψ|² 0.82); VB top interior N p |
| Interior gap (interior-like VB top → CB bottom) | 0.998 eV (Γ → Γ) | ≈ 1.51 eV (Γ: VB → first interior CB) |

- **The H caps remove the Sc-s edge state.** No band near the CB bottom has more than 0.29 edge weight; the CB is
  Sc d on the interior.
- **They add a valence edge band instead.** It disperses 1.24 eV along the edge, from −1.24 eV at Γ to the VBM at X̄
  (doubly degenerate there, like every band at X̄ in this cell), and sits **0.70 eV above the interior VB top**
  (N p, −0.696 eV at Γ). In the [100] H ribbon the Sc–H band is flatter (0.20 eV) and 0.27 eV above the interior.
- So each [110] chemistry has one clean band edge: bare → interior VB, edge CB; H → interior CB, edge VB.
- States at +1.06 / +1.54 / +2.12 eV (Γ) with |ψ|² 0.82–0.84 and Sc s 0.29 are spread over the ribbon faces, not the
  edges (edge 0.02–0.13): vacuum-reaching sheet-like states, as in the bare ribbons.
- **Band survey at 20 Å:** N = 201 / 208 / 217 at L / L + 0.5 / L + 1.0 (25.1 / 26.0 / 27.1 per f.u.; bare [110] at
  15 Å: 21.8 per f.u.). The cost model's 168–200 was low; `cost_a4.py` now uses 201–217 for this ribbon
  (HSE-W90 + PBE-W90 233–651 core-h).
- nbnd 316 was more than needed: band 228 already reaches VBM + 11.68 eV, above L + 1.0 (≈ 10–15 core-h extra).
- Criterion 2 not run (H ribbon, as registered); criterion 5 does not apply (built at 20 Å; whether 20 Å is enough
  is the open 25 Å question).

**Decision (2026-11-13):** the [110] H ribbon gets HSE-W90 in A4b, which now covers four calibration ribbons
([100] n = 3 bare and H, [110] n = 4 bare and H).

## 25 Å check (job 142949, 2026-11-13; cn004, 44 min, 47 core-h; `parse_a4a.py stage2` → `c5_vacuum25`)
NM SCF at 25 Å (`make_a4a.py vac25`) on the two ribbons whose edge pair failed 15 → 20 Å. Same 5 meV gap limit.

| Ribbon | Gap 20 Å → 25 Å (eV) | Δ 20 → 25 | Δ 15 → 20 | Γ shifts, median-aligned (meV) | Verdict |
|---|---|---|---|---|---|
| [100] n = 3 H | 0.7849 → 0.7904 | **5.5 meV** | 17.7 | Sc–H pair (bands 99, 100) −5.6 / −5.4; CB +0.15 | **fails** |
| [110] n = 4 bare | 0.3961 → 0.3911 | **5.0 meV** | 12.3 | Sc-s edge band (CB, 65) −4.85; VB ≤ 0.35 | **on the limit** |

- [110]: QE prints eigenvalues to 0.1 meV, so 5.0 meV is the limit within the data's precision. The script prints
  "fail" only through float rounding (5.0000000000001). Treated as borderline, not as a pass.
- Again only the edge band moves; interior bands ≤ 0.4 meV. ΔE_tot −5.8 / +2.5 meV.
- **Decay:** the 20 → 25 Å step is 0.31 ([100] H) and 0.41 ([110]) of the 15 → 20 step. The earlier 1/L^p guess
  (README above: 6–11 / 4–7 meV) was close: a power law through the three points gives p ≈ 3.5 and 2.5.
- **Residual against infinite vacuum** (two models, three points each; neither is tested):

  | | at 20 Å | at 25 Å | next step 25 → 30 Å |
  |---|---|---|---|
  | [100] H, exponential / power law | 8.0 / 10.1 meV | 2.5 / 4.6 meV | 1.7 / 2.2 meV |
  | [110] bare, exponential / power law | 8.4 / 11.7 meV | 3.4 / 6.7 meV | 2.0 / 2.5 meV |

- So 25 Å would pass the next 5 meV step, but would still leave 3–7 meV on the edge band. 20 Å leaves 8–12 meV.
- **Options (as listed before the run, for a decision):**
  1. accept 20 Å and carry the 8–12 meV edge-band residual into the error budget;
  2. 25 Å for the edge-state ribbons ([100] H, [110] bare, and by the same rule [110] H and the [110] n = 8 hold-out),
     which means redoing A4b-1 for them; project 5.3–14.3k instead of 4.8–12.2k core-h;
  3. keep 20 Å and apply a stated edge on-site correction extrapolated from 15/20/25 Å.

**Decision (2026-11-13): 20 Å stays (option 1).** The edge-band residual (8–12 meV at 20 Å, extrapolated) is carried
in the error budget.
