# T1 corrections: Fröhlich ballisticity and the tunnelling guard (laptop, 2026-11-16)

Settings: `fits/PREREG_corr.md`, frozen at 23:10 on 16 Nov before any device number with them (md5
8aa25f25b7d9313a08dbecf49329e380; no amendment). Code md5s were recorded at 23:11, before the first device run
(`fits/PREREG_corr.code.md5`).
- **Code:** `scnfet/frohlich.py`, `scnfet/guard.py` (tests in `scnfet/tests/test_corr.py`; 68 tests pass) and
  `scripts/t1_corr.py`. Figure and table: `scripts/t1_figs.py corr` → `figures/Fig7b_Ion_W_corrected.png`,
  `summary.md`.
- **Outputs:** one JSON per width (`<family>_n<n>.json`), `dft.json` (the four calibration ribbons), `run.log`.
- **Resources:** one process, peak 0.29 GB, 6 min 43 s. The [100] n = 3 file is from the smoke test, run with the
  same frozen code.
- **Checks passed for every ribbon, carrier and geometry:**
  - the root-found OFF pin reproduces the stored ballistic I_ON within 0.5 %;
  - the re-pin hits I_OFF / R within 1 %.

## Results
### Ballisticity B (L = 15 nm)
B uses the 2D Fröhlich model in vacuum (main); the bulk 3D rate gives the other end.

| | n-FET | p-FET |
|---|---|---|
| [100] bare | m 0.19–0.27; B 0.55–0.63 (3D 0.78–0.84) | m 0.18–0.19; B 0.63–0.64 (3D 0.84) |
| [110] bare | m 0.89–0.98, the Sc-s edge band; **B 0.29–0.31** (3D 0.50–0.53); α 1.03–1.09 | m 0.22–0.39; B 0.47–0.60 (3D 0.72–0.82). n = 4: flat VB (m 38), B 0.04 |
| [110] H | m 0.19–0.29; B 0.54–0.64 | m 0.85–0.86, the Sc–H edge band; **B 0.31–0.32** (3D 0.54); α 1.01–1.02 |

**The heavy edge bands lose about half again as much current to Fröhlich backscattering as the light interior
bands.** Their coupling constant α is just above 1, so the registered flag applies: the perturbative B is not
quantitative there. The direction is robust: heavier means slower and more strongly coupled.

### IRDS at EOT 0.6 after the B correction
Ribbons passing / total. Columns: parallel plate / single gate / GAA. "Ballistic" is before B, from the T1 and
Laplace READMEs.

| | n-FET, ballistic | n-FET, × B | p-FET, ballistic | p-FET, × B |
|---|---|---|---|---|
| [100] bare | 12 / 12 / 12 | **10** / 12 / 12 | 12 / 12 / 12 | **9** / 12 / 12 |
| [110] bare | 0 / 17 / 17 | 0 / **4** / **8** | 15 / 17 / 17 | **0** / 16 / 16 |
| [110] H | 17 / 17 / 17 | **15** / 17 / 17 | 0 / 17 / 17 | 0 / **6** / **11** |

- **Parallel plate × B:**
  - **[100] still meets IRDS at most widths.** The n-FET misses at 1.1 and 1.5 nm (1.52, 1.79). The p-FET
    sits on the line at 5.0–5.9 nm (1.96–1.98 against 1.979).
  - **[110] bare fails both carriers everywhere.** n: 0.40–0.47. p: 0.06–1.45, with the narrowest at 0.06 on a flat
    band.
  - **[110] H keeps the n-FET from 2.0 nm (1.49–2.37) and loses the p-FET everywhere (0.45–0.57).**
- **Laplace single gate × B: the design rule comes back.**
  - [100] passes for both carriers at every width (n 3.75–9.27, p 4.57–11.38).
  - [110] bare n passes only at 1.07–2.0 nm (4 of 17) and falls to 1.16 at 6 nm.
  - [110] H p passes only at 1.4–2.9 nm (6 of 17) and falls to 1.35 at 6.4 nm.
  - The light-band carriers ([110] bare p, [110] H n) pass at every width except [110] bare p at n = 4 (flat band).
- **GAA × B** is in between: [110] bare n passes at 8 of 17 widths, [110] H p at 11 of 17.
- **With the bulk 3D B** (the optimistic end), the Laplace geometries pass everywhere except the flat-band
  point. Under parallel plate, [110] bare and the [110] H p-FET still fail everywhere.
- **Calibration DFT ribbons**, parallel plate × B, then single gate × B (n / p):
  - `r100_3`: 1.52 / 2.09, then 9.23 / 11.34;
  - `r110_4`: 0.39 / 0.08, then 3.12 / 0.56;
  - `r100_3h`: 1.58 / 0.15, then 7.97 / 0.94 (Sc–H VB m 2.2, α 1.6);
  - `r110_4h`: 1.37 / 0.45, then 6.90 / 3.08.

  These sit on the building-block curves.

### Source–drain tunnelling (OFF state, λ_s = 2 nm)
- **R = 1.12–1.38** for the light bands ([100], [110] H n, [110] bare p); 1.00–1.05 for the heavy edge bands.
  **No ribbon reaches the R > 1.5 flag.**
- **Cost:** ΔV_T ≤ 8.3 mV, and re-pinning costs ≤ 2.83 % of I_ON (largest: [110] H n = 14, n-FET, GAA; corrected 2026-11-19 from "≤ 2 %"). The "full" column in `summary.md` is B × re-pinned
  I_ON.
- **At λ_s = 1 nm R is 0.88–0.99:** reflection above a sharp barrier outweighs tunnelling. At 3 nm it is up to
  1.93 for [100] n.

### Band-to-band tunnelling
- **Elastic path through the device:** none. The smallest HSE gap is 0.877 eV, against qV_DD = 0.7 eV, a margin of
  ≥ 0.18 eV.
- **Drain-side window (GIDL-type) upper bound, λ_s = 2 nm:**
  - **[100]: zero at every width.** The window is closed by 0.48–0.77 eV.
  - **[110] bare: flagged at every width from 1.38 nm** (n ≥ 5). The bound is up to 270 × I_OFF. Only n = 4
    (E_g 1.03 eV) is clear.
  - **[110] H: flagged from 2.63 nm** (n ≥ 8; E_g ≤ 1.06 eV), up to 190 × I_OFF.
  - **What the bound means:** it counts every generated carrier as leaving the channel. A floating ribbon has no
    body contact, so the real leakage is lower, but its floating-body behaviour is outside the model. The flag
    says the OFF state of these ribbons is not protected by their gap at V_DD 0.7 V.
- **Registered gap criterion (E_g ≥ 0.9 eV):** [110] bare fails from 2.9 nm (n ≥ 10). [110] H (≥ 0.902 eV) and
  [100] pass.
- **Thermal-tail case:** in some n-FETs the window is closed by a few kT but the bound still flags, because the
  degenerate drain (δ_S up to 0.46 eV) leaves thermal vacancies. An example is [110] H n = 8: window −0.07 eV,
  bound 1.7 × I_OFF. The registered decision quantity is the bound, so these stay flagged.

## What it means for the paper
1. **The orientation/termination design rule survives every correction, and the Fröhlich correction strengthens
   it.** The heavy edge bands that carry it ([110] bare CB, [110] H VB) have B ≈ 0.3 against ≈ 0.6 for the light
   bands.
   - With parallel plate × B, the edge-band carriers fail as in T1. In addition, [110] bare p is lost everywhere,
     [100] n misses at ≤ 1.5 nm, and [100] p sits on the line at 5–6 nm.
   - With the Laplace C_G, the edge-band carriers fail again above about 2–3 nm.
2. **[100] bare is the only family with a clean OFF state.** It has no interband window, no SDT flag, and passes
   IRDS for both carriers at almost every width under every C_G.
3. **[110] bare is weak twice.** Its n-FET is the weakest after B, and above 1.4 nm its OFF state has a drain-side
   interband window. That second point adds to the T1 gap flag rather than replacing it.
4. **Limits to state:**
   - α > 1 for the heavy edge bands;
   - the sheet rate without the ribbon form factor below 2 nm;
   - remote HfO2 phonons are left out;
   - the B factor is the signed-off λ/(λ + L). Fixed top-of-barrier charge would give λ/(λ + 2L), which is lower
     still;
   - the GIDL bound assumes carrier removal.
