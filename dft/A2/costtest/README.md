# HSE cost tests on the A2 sheet (before any ribbon HSE)

Approved 2026-11-11. Written 2026-11-11, **before any output of these tests exists**. They set the
A4b/A5 HSE cost range (`dft/A4a/cost_a4.py`) and the ribbon Wannier outer window. Nothing here enters a fit.

## T1: does smearing remove the empty-band EXX cost?
Our nine A1/A2 HSE runs with `occupations='fixed'` scale as nbnd² (empty bands included): within A2, the per-pair
EXX time rose ×2.45 from nbnd 24 to 56. QE's EXX loops over occupied bands only when it can tell they are empty;
whether smearing changes that is untested.

- **Runs** (`make_costtest.py t1`), A2 sheet, the `hse_k8_q4` settings: HSE06, k 8×8×1 IBZ (15 k), q 4×4×1,
  80/320 Ry, `ecutfock` 160, conv_thr 1e-9, 64 ranks, `-nk 8`, all in one job on one node:

  | Tag | Occupations | nbnd |
  |---|---|---|
  | `ct_fix24` | fixed (same-node control of A2 `hse_k8_q4`) | 24 |
  | `ct_fix40` | fixed | 40 |
  | `ct_sm24` | gaussian 0.005 Ry | 24 |
  | `ct_sm40` | gaussian 0.005 Ry | 40 |

- **Measure:** t = `vexx` wall per call (one call = one k of the pool, one outer step). Ratios r_fix = t(fix40)/t(fix24),
  r_sm = t(sm40)/t(sm24). Quadratic predicts (40/24)² = 2.78, linear (nbnd·nocc) 40/24 = 1.67; the boundary is
  their geometric mean, **2.15**.
- **Decision:**
  - r_fix ≥ 2.15 confirms the nbnd² law inside one system on one node.
  - **r_sm < 2.15 and the smearing eigenvalues equal the fixed ones** (all bands ≤ 24, all 15 k, max |ΔE| ≤ 1 meV,
    after aligning the VBM) → smearing removes the empty-band cost. Then A4b/A5 are staged with gaussian 0.005 Ry
    (it also survives a metallic edge state, where fixed occupations fail) and `cost_a4.py` uses nbnd·nocc for them.
    This is a settings change from "fixed occupations", so it goes for a decision with A4b, not adopted silently.
  - Otherwise fixed occupations and the nbnd² law stay.

## T2: how low can the Wannier outer window go?
Every band under the outer window must be in the HSE run, and HSE cost goes as nbnd² (T1 may change that). A2 used
outer = VBM + 13.5 eV (HSE) / + 11.8 eV (PBE twin), the bulk A1 values.

- **Runs** (`make_costtest.py t2`): wannier90.x re-runs on the A2 `.amn/.mmn/.eig` (16×16×1, 16 projection-only WFs,
  frozen window unchanged at 2.861 eV), with the outer window lowered to what an HSE run with **nbnd = N** would
  allow: outer = min(A2 value, min_k E_N − 0.05 eV), the A2 cap rule. N = 34, 36, 40, 44 (17, 18, 20, 22 bands per
  f.u.; N < 34 leaves fewer than 16 bands at some k). PBE N = 44 equals the A2 baseline and is not rerun.

  | N | HSE outer − VBM (eV) | PBE outer − VBM (eV) |
  |---|---|---|
  | 34 | 8.82 | 8.00 |
  | 36 | 9.74 | 8.44 |
  | 40 | 11.05 | 10.23 |
  | 44 | 12.71 | baseline |
  | A2 baseline | 13.50 (≡ N 48: the A2 run carried 8 bands it did not need) | 11.80 (≡ N 44) |

- **Measure, per variant:**
  1. PBE twin against `pbe_bands` (A2 criterion 4): RMS and max inside the frozen window; max for the VB
     (E ≤ VBM + 0.01 eV) and the CB up to CBM + 1 eV.
  2. HSE against the HSE 16×16 grid eigenvalues inside the frozen window (A2 criterion 5).
  3. HSE against the baseline HSE Wannier bands on the same path: max for the VB and the CB up to CBM + 1 eV.
  4. WF centres (offset from their atoms; z offset of the Sc WFs), spreads, Ω_total.
- **Acceptable N:** (1) VB and CB ≤ CBM + 1 eV within **2 meV** and frozen-window RMS ≤ **3 meV**; (2) ≤ **1 meV**;
  (3) ≤ **5 meV**. The **smallest acceptable outer level** (HSE, relative to the VBM) becomes the ribbon rule for
  A4b/A5; each ribbon's nbnd then follows from its own PBE band survey.
- **Descriptive (no threshold):** whether a lower window pulls the Sc z² WFs back from the vacuum (A2: 0.8–0.94 Å).
  This test only moves the outer window; it does not choose the surface Wannier recipe.

## Run
```bash
python3 make_costtest.py          # writes inputs/, w90/<seed>_nb<N>/, windows.json (local)
# copy costtest/ to /path/to/scn_fet/A2/costtest/ (the W90 runs link ../../../w90/<seed>/)
sbatch run_costtest.sh            # T2 (concurrent serial W90 runs), then T1
# copy outputs back, then:  PYTHONPATH=../../.. ~/.venvs/scnfet/bin/python -P parse_costtest.py
```

**Cost:** T1 ≈ 8–10 min of pw.x on one node, T2 ≈ 4 min, so ≈ 12–16 core-h. That is 2–3× the ≈ 5 core-h quoted
on 2026-11-11, because a same-node control (`ct_fix24`) and a smearing pair were added so each ratio comes from one node.

## Results (2026-11-12; job 142071, cn007, 17 min, 18 core-h; `costtest.json`)
- **T1:** r_fix 2.90, r_sm 2.91: both quadratic. **Fixed occupations stay; the nbnd² law stands.** The eigenvalue
  check gave 1.8 meV, all in band 24 (the top Davidson band) at one k; bands 1–23 agree to 0.1 meV. `ct_fix24`
  ran its EXX ×1.68 slower than the identical A2 `hse_k8_q4` on another node, so C = 2.8–5.9 ns.
- **T2:** N 34 and 36 fail (Wannier functions spread out, ~1 eV errors at the band edges); **N 40 passes**
  (PBE-twin RMS 1.46 meV, VB 1.07, CB 1.03; HSE 1.34 meV from the baseline). **Ribbon rule: outer = VBM + 11.05 eV
  (HSE).** The lower window does not pull the Sc z² WFs back (0.82 Å at N 40).
- Full results: `costtest.json` (T1, T2).
