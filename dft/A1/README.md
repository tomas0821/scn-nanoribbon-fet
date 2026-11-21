# A1: bulk ScN HSE06 + Wannier90, and the PBE twin

Plan §4, stage A1. The criteria below were written down on 2026-11-11, **before any A1 output existed**.

## Settings
- **Structure:** rocksalt ScN, a = **4.501 Å** (experiment; decision 2026-11-11). PBE a₀ from A0 is 4.514 Å (NC).
- **Pseudopotentials:** the all-NC set, Sc SG15 (`Sc_ONCV_PBE-1.0`, 3s3p semicore) + N Dojo, at **80/320 Ry** (A0).
- **HSE06:** `input_dft='hse'`, q 4×4×4, `ecutfock` 160, fixed occupations, nbnd 30, `diago_full_acc`.
- **Wannier90 (stage A1-2):** 8 projection-only WFs (`num_iter = 0`), Sc d + N p, bands 1–5 excluded (Sc 3s3p,
  N 2s). Disentanglement windows are set from the A1-1 eigenvalues. Grid 8×8×8.

## Stage A1-1 (`run_a1.sh`, one 64-rank job, 8 tags)
| Tag | What | Why |
|---|---|---|
| `pbe_scf` | PBE SCF 12³, smearing | reference E_F, density |
| `pbe_scf_bands` | PBE bands L–Γ–X–W–K–Γ | Wannier-quality reference for the PBE twin |
| `pbe_nscf` | PBE NSCF, full 8³ list, nosym | pw2wannier90 input, PBE twin |
| `hse_ibz` | HSE06 SCF, 8³ IBZ | route O |
| `hse_ibz_open` | `open_grid.x` on `hse_ibz` | route O: unfolded to the full grid |
| `hse_nosym` | HSE06 SCF, full 8³ list, nosym | route N (reference) |
| `hse_q8` | HSE06, IBZ, q 8³, nbnd 12 | q convergence |
| `hse_fock320` | HSE06, IBZ, `ecutfock` 320, nbnd 12 | `ecutfock` convergence |

**Cost estimate:** ≈ 35–100 core-h, ≈ 0.5–1.5 h on 64 ranks (kill-panel F2: nosym 8³ q4 nbnd 24 = 13–36 core-h,
×1.25 for nbnd 30; the IBZ runs are a few core-h each; q8 is the largest unknown). Plan: 15–100.

## Pre-registered criteria
1. **Route O vs route N:** HSE eigenvalues of bands 1–20 at all 512 grid points agree within **2 meV**. If so,
   `open_grid.x` (×1.6 cheaper) is used for the A4b ribbons (kill-panel F1); otherwise the nosym route.
2. **q convergence:** the HSE Γ–X, X–X and Γ–Γ gaps change by ≤ **10 meV** from q 4³ to q 8³. If so, q 4³ (and
   the matching sheet q grid) stands; otherwise use q 8³ for bulk and re-cost A2/A3.
3. **ecutfock:** the same three gaps change by ≤ **10 meV** from 160 to 320 Ry. If so, 160 stands.
4. **PBE-twin Wannier quality:** the W90-interpolated PBE bands against `pbe_scf_bands` along the path, inside the
   frozen window: max error ≤ **10 meV**, RMS ≤ **3 meV**. WF centres within 0.05 Å of their atoms.
5. **HSE Wannier:** routes O and N give the same H(R): max |ΔH| ≤ **1 meV**. Interpolated bands reproduce the
   HSE grid eigenvalues inside the frozen window within 1 meV.

**Physics flags** (reported, never tuned to):
- HSE Γ–X ≈ 0.85–0.92 eV (Kumagai 0.88 with semicore at 4.498 Å; Deng 0.92 at 4.50 Å; Mu 0.79 has the semicore
  frozen);
- X–X ≈ 1.9–2.05 eV, Γ–Γ ≈ 3.6–3.8 eV;
- from the W90 bands, m_l ≈ 1.4, m_t ≈ 0.17–0.24;
- PBE Γ–X ≈ −0.03 eV (A0: −0.033 at 4.50 Å).

## Run
```bash
python3 make_a1.py                      # local: writes inputs/ and tags.txt
# copy dft/A1/ to /path/to/scn_fet/A1/, then on the cluster:
sbatch run_a1.sh                        # restartable: finished tags are skipped; tmp_ dirs are kept for stage 2
```

## Stage A1-2 (`run_a1_w90.sh`, after A1-1; written 2026-11-11 after reading the A1-1 eigenvalues)
A1-1 outcome that sets this stage (`parse_a1.py` → `a1_stage1.json`):
- criterion 1 (routes O vs N) **passed**: 0.0 meV at all 512 k-points;
- criterion 2 (q) **failed**: q 4³ → 8³ moves Γ–X by −12.4 meV and X–X by −15.4 meV (> 10). By the rule above,
  bulk production uses **q 8³**, and A2/A3 are re-costed;
- criterion 3 (`ecutfock`) **passed**: 0.0 meV.

Runs:
- `hse_ibz_q8` (HSE06 IBZ, q 8³, nbnd 30) + `hse_ibz_q8_open` (`open_grid.x`): the production HSE06 states;
- Wannier90 (8 projection-only WFs, Sc d + N p, bands 1–5 excluded, `use_ws_distance = .false.`) for:
  `pbe` (PBE twin, from `pbe_nscf`), `hse_n` (q4 route N), `hse_o` (q4 route O), `hse_q8` (q8 route O, production).
- Windows from the A1-1 survey: frozen ≤ VBM + 6.5 eV (HSE) / + 5.4 eV (PBE), below the first non-model band
  (band 14 minimum +7.02 / +5.88 eV); outer ≤ VBM + 13.5 / + 11.8 eV, below band 15 (+13.9 / +12.2 eV).
- Criteria 4 and 5 above apply: `hse_n` vs `hse_o` for criterion 5, `pbe` against `pbe_scf_bands` for criterion 4.
