# A3: the 3-plane slab hold-out (plan §4 A3; §3 D3 thickness axis, LK7; claims gate (b))

**Written 2026-11-15, before any A3 input or output exists.** The plan gives A3 one line ("3-plane slab PBE + HSE + W90
(hold-out)", ≈ 65–180 core-h) and no slab-specific pass criteria. Everything below is a draft until approved. From then
on it is frozen, and later changes are dated amendments below, never edits.

The slab is a pre-registered hold-out: **never refitted**. It tests the SK model along the thickness axis (2 → 3
planes), as the [110] n = 8 and [100] n = 6 ribbons test it along the width.

## Structure
- **3-plane rocksalt-(001) slab**, primitive in-plane cell, 3 f.u. (Sc3N3, 6 atoms).
  - Each plane alternates Sc and N; each plane swaps species relative to the one below.
  - Built by `scnfet.structures.slab(3, …)`, with `frame` R45 as the sheet.
- **Convention as A2 and the ribbons:** a_p = 3.1093 Å (the A0 PBE sheet scaled by bulk a_exp / a_PBE), starting
  plane spacing 2.1685 Å.
- **Vacuum 20 Å** (c = 24.337 Å): the sheet had 20.3 Å and passed its 5 meV check (0.3 meV).
- **Symmetry:** mirror-symmetric about the middle plane, so it has no dipole and needs no dipole correction (as the
  sheet).

## Stages
**A3-1, PAW PBE relaxation** (as A4a / A5: Sc PAW psl 0.2.3 spn, N Dojo; 90/720 Ry; gaussian 0.005 Ry):
- Internal coordinates at fixed a_p and c. By symmetry only z moves (plane spacing and rumpling). BFGS,
  `forc_conv_thr` 1e-3, `etot_conv_thr` 1e-5; k 12 × 12 × 1 (A0's sheet grid); nspin 1.
- **No magnetism check.** Bulk, the sheet and every ribbon edge (n = 3–8) are non-magnetic, and the slab cuts no
  in-plane bond.
- One small job (`run_a3.sh`, = A5's `run_a5.sh` with A3's directory).

**A3-2, NC PBE + HSE06, as A2** (SG15 Sc, Dojo N, 80/320 Ry), staged from the relaxed geometry, one job:
- **PBE:**
  - SCF 16 × 16 × 1 (A2: m-v 0.005 Ry);
  - `pp.x` bare + Hartree potential → vacuum level;
  - bands Γ–X̄–M̄–Γ, 60 points per segment, with `projwfc.x` on them (the W2 character test);
  - NSCF on the full 16 × 16 × 1 list (nosym) for Wannier90;
  - **vacuum check:** SCF with vacuum 27.5 Å (c + 7.5 Å, as A2's 22.5 → 30), 12 × 12 × 1.
- **HSE06:**
  - IBZ SCF 16 × 16 × 1, q 8 × 8 × 1 (A2's production q), fixed occupations, `ecutfock` 160;
  - `pp.x` (vacuum level) and `projwfc.x` (non-model states, H1 / H3);
  - `open_grid.x`, then `pw2wannier90` (route O, as A1 / A2).
- **Band count 84** in the NSCF and the HSE run: A2's 28 per f.u. (56 for the sheet). The windows criteria check that
  the outer windows are not capped.

**A3-3, Wannier90, the A4b recipe as frozen on 13 Nov plus its two amendments of 14 Nov** (`dft/A4b/README.md`):
- **Set B** = Sc d + N p (24 WFs; Sc 3s 3p + N 2s excluded, 15 bands). **Set S** = B + Sc s on every Sc (27 WFs).
- Projection-only (`num_iter` 0). Outer window VBM + 10.23 eV (PBE) / + 11.05 eV (HSE), the T2 rule.
- **Frozen windows:** W1 (A2's rule); W2 (W1 and 0.1 eV below the first non-model state, s_out ≥ 0.25).
- **W3 on set S, under amendment 1's rule only:** if an HSE non-model state lies inside [VBM − 0.5, CBM + 0.5],
  W3 is used if it passes H1, H2 and covers that window; otherwise W2 is used with a stated limitation. As for the
  hold-outs, this is reported, not a stop.
- **Decision rule as A4b:** W2 → W1 → back to the authors.

## Run-level criteria (judged when the runs are in)
- **Vacuum (A2 criterion 2):** the PBE gap and the M̄ offset change by ≤ 5 meV from 20 to 27.5 Å.
- **P1–P5 on the PBE twin, as frozen** (device window [VBM − 0.5, CBM + 0.5]).
  - The reference is the NC PBE bands on the Γ–X̄–M̄–Γ path.
  - P5's masses are taken at both band edges along Γ–X̄ and Γ–M̄ (the 2D analogue of the ribbons' single direction).
  - P6 is descriptive (WF centres; Sc z² offsets; with S, the Sc s offsets).
- **H1, H2 (≤ 1 meV on the HSE grid inside the frozen window), H3 (flag; amendment 1's rule as above).**
- **H4 (amendment 2, as for the hold-outs):** variants that freeze the same window agree within 5 meV over the
  device window on a fine path (200 points per segment). Otherwise it goes back to the authors.

## Hold-out criteria (claims gate (b); judged at the gate, not here). Proposed, as D7's ribbon numbers without the device metrics
- **PBE-twin SK** (refit to PBE bulk + sheet only, never to the slab) **vs the PBE-Wannier slab:**
  - gap ±50 meV;
  - valley offsets E_c(M̄) − E_c(Γ) and E_c(X̄) − E_c(Γ) ±50 meV;
  - in-plane transport masses at both band edges ±10 %.
- **SK-HSE vs the HSE slab:** gap ±80 meV (D7's HSE single-point limit).
- **Failure → plan LR2:** no refitting; reported as a transferability limit along the thickness axis.
- **Descriptive:**
  - gap and M̄ offset against the sheet (the thickness trend; LR4 context);
  - the position of the vacuum-like state (sheet: HSE CBM + 0.74 eV);
  - surface character of the band edges;
  - the relaxed plane spacings, against the sheet's 2.17 Å and bulk's 2.25 Å.

## Cost
| Stage | core-h |
|---|---|
| A3-1 relaxation | ≈ 3–10 |
| A3-2 PBE (≈ 5) + HSE06 | ≈ 140–260 |
| A3-3 Wannier90 | ≈ 2–10 |
| **A3 total** | **≈ 150–280** |

The HSE figure is the sheet's measured 68 core-h (A2 `hse_ibz_k16`, 64 min, nbnd 56) × (84 / 56)² × 1.08 (c), with
a range of ×0.8–1.5. The plan said ≈ 65–180 and `cost_a4.py` 100–200. The difference is the band count: A2's 28 per
f.u. at 3 f.u.

## Status
- **2026-11-15: frozen on approval.** `README.prereg.md5` holds the md5 of lines 1–84 (everything above
  this section). A3-1 submitted: job 143638 (`relax_slab3`, cn009, 14:01).

## Results: A3-1, the slab relaxation (2026-11-15)
- **Job 143638 (cn009):** 1 min 54 s, ≈ 2 core-h. BFGS converged in 6 steps (7 SCF); final total force 0.0002
  Ry/Bohr.
- **Geometry (descriptive):**
  - The outer planes move out by 0.06–0.08 Å.
  - Spacing from the outer to the middle plane: 2.231 Å (outer Sc → middle N) and 2.245 Å (outer N → middle Sc),
    against the sheet's 2.17 Å and bulk's 2.25 Å.
  - The outer N sits 0.014 Å outside its Sc (rumpling); the middle plane does not move (mirror symmetry kept).
  - Relaxation energy 54 meV per cell.

## A3-2 staged (2026-11-15)
`make_a3.py stage2` → `inputs/{scf,pot,bands,projwfc,nscf,hse,hpot,hproj,hopen,vac}_slab3.in`,
`w90/slab3/slab3_{B,S}.{win,pw2wan}`, `w90/slab3/H/slab3_{H,HS}.{win,pw2wan}`, `tags_slab3.txt`, `tags_hse_slab3.txt`,
`tags_vac.txt`.
- **Positions are the relaxation's.** The settings match A2's (`hse_ibz_k16`: the header differs only in the cell
  form, nat and nbnd 84).
- **Three jobs, all in `scn_fet/A3`:**
  1. **PBE twin:** A4b-1's `run_a4b.sh`, unchanged (`A4B_DIR`, R = `slab3`; the tag names are A4b-1's). Minutes.
  2. **HSE06:** `run_a3_hse.sh`, i.e. `run_a4b2.sh`'s steps with `run_a5_hse.sh`'s convergence guard. It runs on the
     45 IBZ k with 63 ranks in 9 pools, checks the open_grid list for 256 points, and does Wannier90 `-pp` +
     pw2wannier90 for H and HS. ≈ 2–4 h; 12 h limit.
  3. **Vacuum check:** `run_a3.sh vac_slab3`, the relaxed slab centred in 27.5 Å of vacuum (c 31.837 Å),
     12 × 12 × 1, nbnd 36 (A2's `pbe_c30` had 24 for 2 f.u.). Minutes.
- **The jobs don't interfere:** their tmp dirs, outputs and Wannier files are disjoint. The PBE twin's clean-up
  globs do not reach `w90/slab3/H/`.
- **Stub dry run:** all 8 + 6 + 1 tags in order, reruns skip all, pools right, both HSE `.win` files filled with the
  256-point list.
- **Cost ≈ 145–270 core-h**, almost all of it the HSE run, within the frozen estimate.

## A3-2 submitted (2026-11-15, approved)
24 files md5-checked into `scn_fet/A3` (`run_a4b.sh` byte-identical to A4b's). Submitted 14:36, all started at once:
**143662** PBE twin (cn009), **143663** HSE06 (cn012; 63 ranks, 9 pools, nbnd 84), **143664** vacuum check.

## Results: vacuum check (job 143664, 13 s) — passes
PBE, 16 × 16 at 20 Å against 12 × 12 at 27.5 Å; Γ, X̄ and M̄ lie on both grids (`a3_vacuum.json`):

| | 20 Å | 27.5 Å | Change |
|---|---|---|---|
| Gap (direct at Γ) | 0.1201 eV | 0.1190 eV | −1.1 meV |
| M̄ offset E_c(M̄) − E_c(Γ) | 0.6217 eV | 0.6219 eV | +0.2 meV |
| E_c(X̄) − E_c(Γ) | 3.1426 eV | 3.1429 eV | +0.3 meV |

Both registered quantities are inside 5 meV, so 20 Å stays.

**Descriptive (PBE):** the slab gap is 0.120 eV, direct at Γ (sheet 0.404). The M̄ offset is 0.622 eV (sheet 0.876).
Both fall with thickness towards bulk, which is a PBE semimetal (Γ–X −0.034 eV, A0). This fits the folding:
- in the primitive (001) cell, the bulk X valley along the slab normal folds to Γ̄ and the two in-plane X valleys
  to M̄;
- each added plane brings the folded valleys closer together.

For the valley-crossover question (LR4) this is the thickness axis, not the width axis. HSE moves both (A3-2).

## Results: A3-2 PBE twin and the A3-3 windows (job 143662, 8 min 34 s; 2026-11-15)
`parse_a3.py windows` → `a3_windows.json` and `w90/slab3/{W1,W2,W3}/slab3_<v>.win`. `parse_a3.py` is new: the frozen
A4b window rules, in 2D (`parse_a4b.py`'s constants and projwfc parser; A2's vacuum-level rule for a centred slab).
Its `vacuum` mode reproduces `a3_vacuum.json` exactly (that check was first parsed inline).
- **PBE (NC):**
  - gap 0.120 eV, direct at Γ, on the path and on the 16 × 16 grid;
  - E_vac − VBM 3.304 eV (27 z points ≥ 8.8 Å from the atoms; spread 0.03 meV).
- **First non-model state (set B, s_out ≥ 0.25): CBM + 1.459 eV at Γ** (band 27; s_out 0.55, |ψ|² 0.86; Sc d 0.50,
  Sc s 0.42, N s 0.06).
  - It is the same kind of vacuum-reaching state as the sheet's (PBE CBM + 0.885, HSE + 0.743 eV;
    `dft/A2/a2_sheet_nonmodel.json`).
  - The threshold doesn't matter here: the same state at 0.20 and 0.30. With set S (+ Sc s), the first one is
    CBM + 2.28 eV.
  - **Amendment 1 (PBE): not inside the device window.** Whether W3 runs depends on the HSE projections (A3-3 rule).
- **Windows (frozen rules):**

  | Variant | Set, WFs | Frozen top | Bound by | Vacuum margin | Bands frozen ≤ / in outer ≥ |
  |---|---|---|---|---|---|
  | W1 | B, 24 | CBM + 2.684 eV | vacuum | 0.50 eV | 14 / 30 → ok |
  | W2 | B, 24 | CBM + 1.359 eV | non-model state | 1.82 eV | 11 / 30 → ok |
  | W3 | S, 27 | CBM + 2.684 eV | vacuum | 0.50 eV | 14 / 30 → ok (`.win` written; runs only under amendment 1) |

  - Outer window VBM + 10.23 eV, not capped.
  - The band path covers the highest frozen top + 0.1 eV.
  - Both W1 and W2 freeze the whole device window (CBM + 0.5).

## A3-3 PBE variants submitted (2026-11-15)
W1 and W2 run in A5-2's job 2 (**job 143673**, `dft/A5/twin/run_a5_w90.sh … A3:slab3`; `dft/A5/README.md`). The HSE
variants follow the HSE run (143663) as a separate small job.

## Results: A3-3 PBE variants (job 143673, 2026-11-15)
`parse_a3.py criteria` → `a3_criteria.json`: P1–P5 as frozen, in 2D. The reference is the NC PBE bands on
Γ–X̄–M̄–Γ (181 points); the Wannier bands come from `scnfet.w90.hk_frac`. P5 is taken at both band edges along
Γ–X̄ and Γ–M̄. H4 uses 200 points per segment.

| | W1 | W2 | Limit |
|---|---|---|---|
| P1 (vacuum margin, counts, outer) | ok | ok | |
| P2 device-window max error (86 values) | 1.02 meV | 0.92 meV | ≤ 5 |
| P3 frozen RMS (max) | 0.82 (17.7) | 0.23 (3.0) | ≤ 3 |
| P4 gap | −0.04 meV | −0.05 meV | ≤ 3 |
| P5 VB / CB along Γ–X̄ | −0.1 / −0.8 % | −0.1 / −1.2 % | ≤ 3 % |
| P5 VB / CB along Γ–M̄ | +0.2 / −0.5 % | −0.3 / −1.0 % | ≤ 3 % |
| P6 outer-plane Sc z² offsets | 0.91 Å | 0.19 Å | descriptive |
| P6 other host WFs, max offset | 0.05 Å | 0.09 Å | descriptive |

- **Decision rule: W2.** H4 W1 vs W2: 0.57 meV over the whole device window, ≤ 5 meV.
- **Masses (PBE, DFT path), m_e:**
  - VB (degenerate at Γ, 0.00 meV): upper branch −0.355 along Γ–X̄ and −0.189 along Γ–M̄; second branch −0.116
    and −0.187;
  - CB: 0.099 along Γ–X̄, 0.108 along Γ–M̄.
- **Valleys (PBE):** E_c(M̄) − E_c(Γ) 0.6216 eV (W2 0.6216).
  - **E_c(X̄) − E_c(Γ) = 3.143 eV in DFT, but 3.631 in W2 and 3.150 in W1.** That X̄ state lies above both frozen
    tops (W1 CBM + 2.68, W2 + 1.36). It is a mixed state with s_out 0.33 (Sc d 0.78, Sc s 0.15, N s 0.07), non-model
    by the frozen rule.
  - **For the claims gate (not changed here):** hold-out criterion "E_c(X̄) − E_c(Γ) ±50 meV vs the PBE-Wannier
    slab" has a reference that W2 does not reproduce. The offset is 3.1 eV above the CBM, far outside any device
    window. It needs a decision at the gate, as a dated amendment if anything changes (e.g. the DFT value as the
    reference for that one offset).
- **P6, the centres:** the reported Ω (1,205 Å² in both) is an artifact of the centred slab.
  - The 8 WFs of the middle plane (Sc d ×5, N p ×3) sit at exactly z = c/2. With one k along z, their
    Berry-phase centres fall on the branch cut, so Wannier90 reports z ≈ 0 and a spread of ≈ (c/2)² = 148 Å².
  - Off the cut, the 16 WFs sum to 23.6 Å² (W1) and 18.9 Å² (W2), max 3.2 / 1.2 Å².
  - The Hamiltonian does not use the centres (`use_ws_distance = .false.`, kz = 0 only). `parse_a3.centres` flags
    these WFs and keeps only their in-plane offset.

## Results: A3-2 HSE06 (job 143663, 3 h 17 min, 210 core-h; 2026-11-15) and the HSE windows
`parse_a3.py hwindows` → `a3_hwindows.json`, `w90/slab3/H/{W1,W2}/slab3_H<v>.win`.
- **Converged:** 5 EXX outer steps, dexx 1.3e-10 Ry (`!!`). `open_grid` and pw2wannier90 for sets B and S:
  `JOB DONE`.
- **HSE gap 1.1055 eV, direct at Γ** (PBE 0.120; sheet HSE 1.415).
- **Valleys:** M̄ offset 0.676 eV (PBE 0.622; sheet HSE 0.99); E_c(X̄) − E_c(Γ) 3.317 eV.
- **E_vac − VBM 3.860 eV** (PBE 3.304).
- **First non-model state (set B): CBM + 1.360 eV at Γ**, the vacuum-reaching state (Sc s 0.46, Sc d 0.45, N s 0.07;
  s_out 0.63). The same at 0.20 and 0.30.
  - PBE had it at + 1.459; the sheet's HSE state is at + 0.743. A third plane pushes it up.
  - **H3 / amendment 1: outside the device window, so no W3.**
- **H1: both variants ok.**

  | Variant | Frozen top | Bound by | Vacuum margin | Bands frozen ≤ / in outer ≥ |
  |---|---|---|---|---|
  | HW1 | CBM + 2.254 eV | vacuum | 0.50 eV | 14 / 31 |
  | HW2 | CBM + 1.260 eV | non-model state | 1.49 eV | 11 / 31 |

  - Outer VBM + 11.05 eV, not capped (top band minimum VBM + 19.56).
  - Both freeze the whole device window.
  - The two `.win` files differ from `slab3_H.win` only in the two window lines.

## A3-3 HSE variants staged (2026-11-16)
- **Job script:** `dft/A4b/run_a4b2_w90.sh`, unchanged (md5 equal on the cluster), with `A4B_DIR` pointed at
  `scn_fet/A3`.
  - It links `slab3_H.{amn,mmn,eig}` into `H/W1`, `H/W2` and runs both side by side.
  - Its set-S step is skipped (`slab3_HS.amn` and its `JOB DONE` are there); W3 is skipped (no `.win`).
- **Stub dry run:** both variants started with the right links; a rerun skips them.
- **To upload:** the two `.win` files. Then from `scn_fet/A3`: `A4B_DIR=$PWD sbatch -J a3hw ../A4b/run_a4b2_w90.sh
  slab3`.
- **Cost:** ≈ 10–20 min on one node, ≈ 10–20 core-h.
- **After it:** H2 (≤ 1 meV on the HSE grid inside the frozen window), H4 between HW1 and HW2, decision as
  registered (W2 → W1 → back to the authors).

## A3-3 HSE variants submitted (2026-11-16, approved)
The two `.win` files were md5-checked into `scn_fet/A3`. **Job 143764** on cn003, started 06:58:24, from `scn_fet/A3`:
`A4B_DIR=$PWD sbatch -J a3hw ../A4b/run_a4b2_w90.sh slab3`. Both variants started with the `slab3_H` links and were in
the disentanglement (≈ 0.09 s per iteration). `parse_a3.py hcriteria` (H2, H4, decision) is written for when it ends.

## Results: A3-3 HSE variants (job 143764, 2 min 14 s, 2.4 core-h; 2026-11-16)
Both runs "All done". Disentanglement converged in 1,201 (HW1) / 1,301 (HW2) iterations. Copied back as for the PBE
variants. `parse_a3.py hcriteria` → `a3_hcriteria.json`.

| | HW1 | HW2 | Limit |
|---|---|---|---|
| H1 | ok | ok | |
| H2, full open_grid (256 k), the set Wannier90 froze | 0.042 meV (2,603 values) | 0.032 meV (2,372) | ≤ 1 |
| H2, IBZ stdout eigenvalues (cross-check) | 2.643 meV (0.078 without the state at the top) | 0.074 meV | |
| Gap (fine path) | 1.1054 eV | 1.1055 eV | grid 1.1055 |
| E_c(M̄) − E_c(Γ) | 0.6756 eV | 0.6756 eV | grid 0.6756 |
| E_c(X̄) − E_c(Γ) | 3.3235 eV | 3.7756 eV | grid 3.3173 |
| Outer-plane Sc z² offset | 0.98 Å | 0.31 Å | descriptive |
| Off-cut Ω (16 WFs) | 26.1 Å² | 21.6 Å² | descriptive |

- **Decision: HW2** (W2 → W1 → the authors). **H4 HW1 vs HW2: 0.70 meV** over the whole device window (≤ 5).
- **How H2 is scored:** on the full open_grid, with the `.eig` eigenvalues (12 digits) and `dis_froz_max` as written in
  the `.win`, i.e. the states Wannier90 actually froze.
  - The first parse used QE's stdout IBZ eigenvalues (4 decimals) and the unrounded frozen top. That put one HW1 state
    at 4.2163 eV inside the window. The `.eig` file has it at 4.216311667 eV, 0.012 meV above the `.win`'s 4.2163, so
    Wannier90 did not freeze it, and its 2.64 meV is a disentangled-state error.
  - No decision depends on this: HW2 passes on either count.
  - A4b's `hvariant` (H2 ≤ 0.32 meV, all pass) uses the stdout count. That is only an issue for a state within
    0.05 meV of a frozen top.
- **Masses (HSE, fine path, descriptive), m_e:**
  - CB 0.170 (Γ–X̄) / 0.170 (Γ–M̄); PBE 0.099 / 0.108.
  - VB (degenerate at Γ): −0.342 / −0.184 along Γ–X̄ and −0.374 / −0.177 along Γ–M̄.
  - These come from the parabolic fit within 30 meV of Γ on 200 points per segment. The PBE numbers used the
    61-point DFT path, so the two are not fitted alike.
- **E_c(X̄) − E_c(Γ):** HW2 misses the X̄ state by 0.46 eV, as in PBE (above its frozen top; non-model). It is the same
  claims-gate flag as for the PBE slab.
- The HW2 Hamiltonian `w90/slab3/H/W2/slab3_HW2_hr.dat` is the HSE slab reference for gate (b)'s ±80 meV gap
  criterion (the gap needs only the grid value, 1.1055 eV).
