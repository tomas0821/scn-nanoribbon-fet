# Claims gate (e): the ToB code on the HfS2 zigzag benchmark (laptop, 2026-11-16)

- **Benchmark:** Matić & Poljak, Electronics 13, 1048 (2024), zigzag HfS2 nanoribbons.
- **Conventions and targets:** pre-registered and signed off on 11 Nov (restated in the `scnfet/tob.py`
  docstring), before any benchmark run. Reported, not tuned.
- **Inputs:** QUDEN_MATERIALS 2025_01 `HfS2XNR` Wannier Hamiltonians, N = 2–6, 8 and 10 (Poljak & Matić, CC BY-NC-SA
  4.0; `benchmark/quden/`, md5s in `MD5SUMS`), plus QUDEN's own parameters for the set (`material_selection.m`).
- **Script:** `scripts/gate_e_hfs2_benchmark.py`, writing `gate_e.json`, `run.log` and `figures/Fig9_HfS2_benchmark.png`.
- **Resources:** 2 threads and 2.5 GB caps; peak RSS 0.49 GB; 18 min.

## Width pre-check (registered branch)
- Neither QUDEN package contains structure files, so our-rule W (outermost nuclei) cannot be computed.
- QUDEN's ribbon widths, b = 13.61 … 63.97 Å, are exactly Matić's quoted 1.36 + 0.629 (N − 2) nm.
- **By the pre-registered branch, the benchmark is judged on Matić's W.** That validates the code, not the width
  definition. There is no second HfS2 set to report.

## Input checks
- **Wannier count:** 22N + 6 (Hf d5 and 2 × S p3 per formula unit, plus 2 edge S with p3).
- **Valence count:** QUDEN's CB index − 1 = 3 (4N + 2), the S p bands.
- **Channel atoms for the doping:** 6N + 2 (no H in these ribbons). The S/D doping comes out at
  2.67–2.83 × 10¹³ cm⁻², the paper's "≈ 2.7 × 10¹³".
- **Band gaps:** our 2000-k gaps match QUDEN's listed gaps within 0.0–1.6 meV.
- **k-convergence:** I_ON changes by ≤ 8 × 10⁻⁶ from 2000 to 4000 k.
- **SS:** 59.5 mV/dec, as ideal electrostatics gives.

## Result: **all five registered points are inside ±10 %** (EOT 1 nm, VDD 0.7 V, I_OFF 10 nA/μm)
| N | W (nm) | n-FET I_ON (mA/μm) | paper | Δ | p-FET I_ON | paper | Δ |
|---|---|---|---|---|---|---|---|
| 2 | 1.36 | 1.696 | 1.66 | + 2.2 % | 1.223 | 1.18 | + 3.7 % |
| 3 | 2.00 | 1.738 | 1.74 | − 0.1 % | 1.189 | | |
| 4 | 2.63 | 1.537 | | | 1.212 | | |
| 5 | 3.25 | 1.357 | | | 1.242 | | |
| 6 | 3.88 | 1.198 | | | 1.276 | | |
| 8 | 5.14 | 1.069 | | | 1.345 | | |
| 10 | 6.40 | 1.048 | 1.14 | − 8.1 % | 1.408 | 1.50 | − 6.1 % |

- **Descriptive:** the paper's "pFET … 1.12 at narrow widths" against our minimum over N = 2–4, 1.189: + 6.2 %.
- **The trends match the paper.** n peaks near 2 nm and falls with width; p rises with width; the two cross near
  3.6 nm.
- **At 6.4 nm both of ours sit 6–8 % below the paper.** That is within the band, so there's no LR10 and nothing was
  investigated. Two candidates to name in the paper:
  - QUDEN's own ToB defaults use gate / drain couplings of 0.99 / 0.01 (`main_tob.m`), where our registered
    ideal values are 1 / 0;
  - the precision of reading values off the paper's figure.
- **EOT 0.6 nm** (descriptive; the paper's targets are at EOT 1 nm): n 1.75–2.97, p 2.02–2.49 mA/μm; in `gate_e.json`.

## Verdict
**Gate (e) reported, inside the registered ±10 % band at every registered point** (largest deviation 8.1 %). The ToB
code and the signed-off conventions reproduce the published HfS2 benchmark, so all five claims gates are now run.
