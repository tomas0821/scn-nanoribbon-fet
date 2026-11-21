# A2: the M-ScN sheet, HSE06 + Wannier90, and the PBE twin

Plan §4, stage A2. The criteria below were written down on 2026-11-11, **before any A2 output existed**.
Staged for approval; nothing has been submitted.

## Settings
- **Structure:** M-ScN, P4/nmm two-plane rocksalt-(001) bilayer, 2 f.u. (4 atoms), c = 22.5 Å (A0).
  - The geometry is the A0 PBE (PAW 90/720) 2Dxy vc-relax: a = 3.1179 Å, Sc/N planes 2.174/2.177 Å apart.
  - It is **scaled by the bulk ratio 4.501 / 4.5134** (a_exp / a_PBE): **a = 3.1093 Å**, planes 2.168/2.171 Å,
    in-plane Sc–N 2.199 Å. This is the convention the bulk used (experimental a, the A1 choice), so sheet-vs-bulk
    strain stays as PBE predicts it.
  - Unscaled alternative: the plain PBE geometry (a = 3.1179 Å).
- **Pseudopotentials and cutoffs:** as A1: NC Sc SG15 + N Dojo, 80/320 Ry; HSE06 `ecutfock` 160 (A1 criterion 3).
- **HSE06:** fixed occupations, q 8×8×1. Bulk q 4 → 8 moved the gaps by 12–15 meV (A1 criterion 2 failed), so q 8
  is the starting point and q 4 vs q 8 is re-checked here.
- **Wannier90 (stage A2-2):**
  - 16 projection-only WFs (2 × Sc d + 2 × N p); bands 1–10 excluded (2 × Sc 3s3p, 2 × N 2s).
  - **16×16×1 grid.** Bulk 8³ failed the interpolation criterion and 16³ fixed it; the sheet's in-plane a (3.11 Å)
    is close to the bulk primitive length (3.18 Å).
  - Route O (IBZ + `open_grid.x`, A1 criterion 1: identical to route N).

## Stage A2-1 (`run_a2.sh`, one 64-rank job, 9 tags, restartable)
| Tag | What | Why |
|---|---|---|
| `pbe_scf` | PBE SCF 16×16×1, m-v 0.005 Ry, nbnd 56 | density, E_F |
| `pbe_pot` | `pp.x` bare + Hartree potential → cube | vacuum level: the frozen window must stay below it |
| `pbe_bands` | PBE bands Γ–X̄–M̄–Γ, 60 points per segment | reference for the PBE-twin criterion |
| `pbe_nscf_k16` | PBE NSCF, full 16×16×1 list, nosym | pw2wannier90, PBE twin |
| `pbe_c30` | PBE SCF at c = 30 Å, 12×12×1 | vacuum convergence |
| `hse_k8_q4`, `hse_k8_q8` | HSE06 SCF 8×8×1 IBZ, q 4 / q 8, nbnd 24 | q convergence |
| `hse_ibz_k16` | HSE06 SCF 16×16×1 IBZ (45 k), q 8, nbnd 56 | production |
| `hse_ibz_k16_open` | `open_grid.x` | full 256-k grid for pw2wannier90 |

**Cost estimate** (scaled from the measured A1 runs: k, q, bands, occupied bands and the 9.6× larger cell):

| Step | core-h |
|---|---|
| PBE steps | ≈ 5–10 |
| q check | ≈ 5 + 20 |
| production HSE | ≈ 95 |
| **Total** | **≈ 130 (range 80–250)** |

That is ≈ 2–4 h on 64 ranks, against the plan's 50–100 core-h; the overrun is A1's denser grids. The A0 + A1 + A2
running total is ≈ 290 core-h, against the ≈ 4k anchor.

## Pre-registered criteria
1. **q convergence:** the HSE direct gap at Γ, the M̄ valley offset E_c(M̄) − E_c(Γ) and the X̄ conduction level
   change by ≤ **10 meV** from q 4 to q 8 (k 8×8). If so, q 8 is shown converged against q 4 for the sheet;
   otherwise the shift is reported and q 8 stays, as for bulk.
2. **Vacuum:** the PBE gap and M̄ offset change by ≤ **5 meV** from c = 22.5 to 30 Å.
3. **Windows (A2-2):** the frozen window top is **≥ 0.5 eV below the vacuum level**, taken as the planar-averaged
   potential in the vacuum, and below the first non-model band. The outer window is set from the A2-1 survey, as
   for A1.
4. **PBE-twin Wannier quality:** the W90 bands against `pbe_bands` inside the frozen window: max ≤ **10 meV**,
   RMS ≤ **3 meV**; WF centres within 0.05 Å of their atoms. This is the same test as A1 criterion 4, which bulk
   passed on RMS but failed on the strict max at 16³.
5. **HSE Wannier:** the interpolated bands reproduce the HSE grid eigenvalues inside the frozen window within
   **1 meV**.

**Physics flags** (reported, never tuned to):
- PBE gap ≈ 0.40 eV, direct at Γ (A0: 0.406 at a = 3.118 Å; the 0.27 % compression shifts it slightly);
- HSE06 gap ≈ 1.39 eV, direct (Liang 2021, PBE-D3 geometry, a = 3.11 Å; 1.38 eV in their ref. 27);
- **M̄ valley offset:**
  - PBE at A0's geometry: E_c(M̄) − E_c(Γ) = 1.264 − 0.406 = **0.86 eV**;
  - plan risk LR4 triggers at > 0.6 eV ("valley crossover outside the DFT bracket: report as a model prediction");
  - the HSE value decides LR4 for the sheet.

## Run
```bash
python3 make_a2.py          # local: writes inputs/ and tags.txt (done; dry-run with a stub srun passed)
# copy dft/A2/ to /path/to/scn_fet/A2/, then on the cluster:
sbatch run_a2.sh            # restartable; tmp_pbe_scf and tmp_hse_ibz_k16 are kept for stage 2
```
