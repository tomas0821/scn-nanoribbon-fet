# Post hoc checks, 2026-11-19 (not pre-registered)

Descriptive sensitivities quoted in the manuscript after the pre-submission review. They reuse stored results only
(no new DFT, no new device run) and change no registered number or verdict.

| Output | Script | Quoted in |
|---|---|---|
| `posthoc_b_variants.txt` | `scripts/posthoc_b_variants.py` | main §3.4 (B ranges), SM S4 and Table S2 (pass counts under ε12 = 4, 5 and the charge-pinned form λ0/(λ0 + 2L)) |
| `posthoc_ion_margins.txt` | `scripts/posthoc_ion_margins.py` | Table 3 † and ‡ marks (widths within ±4 % of the IRDS target) |
| `posthoc_hse_holdout_kcheck.txt` | `scripts/posthoc_hse_holdout_kcheck.py` (with `posthoc_qe_parse.py`) | main §3.1, SM S1 (HSE06 building blocks against the A5-3 HSE06 SCF eigenvalues of the hold-outs) |
| `posthoc_pbe_chord.txt` | `scripts/posthoc_pbe_chord.py` | main §3.1, SM S1 (the same Γ→k1 chord in PBE, for comparison) |

Run from the project root: `PYTHONPATH=. python scripts/<script>.py > fits/posthoc_2026-11-19/<script>.txt`.

Key numbers:
- B (L = 15 nm), heavy / light bands: 0.29–0.32 / 0.47–0.64 (registered, ε12 = 1); 0.58–0.61 / 0.80–0.91 (ε12 = 4);
  0.64–0.68 / 0.84–0.93 (ε12 = 5); bulk 3D 0.50–0.54 / 0.72–0.85; charge-pinned 0.17–0.19 / 0.31–0.47.
- HSE06 hold-out check: r110_8 CB ≤ 7.8 meV, VB ≤ 36 meV, chord masses −6.4 / +24.7 % (PBE −3.0 / +7.9 %);
  r100_6 CB ≤ 24 meV, VB ≤ 19 meV, chord masses −0.1 / +1.7 %. Gaps reproduce the registered −9.9 / −24.5 meV.
- 13 cells of Table 3 have a width within ±4 % of the target; the [110] bare p, s = 2 nm cell has 8 of 17.
