# Post hoc check, 2026-11-20 (not pre-registered)

**Nearly-free-electron (NFE) band of bare [100] n = 3.** In HSE06 a vacuum-reaching NFE state sits 0.35 eV above the
conduction-band minimum of this calibration ribbon, inside the device window (dft/A4b/README.md, 14 Nov amendment).
`scripts/posthoc_nfe_share.py` computes the n-FET on-current with HW3 (the DFT device Hamiltonian, which contains the
band) and with HW2 (which leaves it out), at the registered ToB settings and gate capacitances.

| C_G | EOT | HW3 (mA/µm) | HW2 (mA/µm) | NFE share |
|---|---|---|---|---|
| parallel plate | 1.0 | 1.585 | 1.576 | +0.5 % |
| parallel plate | 0.6 | 2.754 | 2.740 | +0.5 % |
| isolated single gate | 1.0 | 15.546 | 15.575 | −0.2 % |
| isolated single gate | 0.6 | 16.762 | 16.807 | −0.3 % |
| isolated GAA | 1.0 | 16.989 | 17.038 | −0.3 % |
| isolated GAA | 0.6 | 18.289 | 18.363 | −0.4 % |

The HW3 parallel-plate value at EOT 0.6 reproduces the stored DFT point (fits/t1/dft.json, 2.7542 mA/µm). The band is
at CBM + 0.351 eV at Γ in HW3 and absent in HW2; the gap is the same in both (1.9643 eV). Quoted in the manuscript's
Limitations and in SM S1. The building-block curves contain no NFE band (neither source carries interior Sc s).
Output: `nfe_share.json`.
