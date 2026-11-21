# T1 with the geometric gate capacitance (Laplace C_G; laptop, 2026-11-16)

Settings: `fits/PREREG_cg.md`, frozen at 20:28 on 16 Nov before any device number with it (md5
473c7583fcaf4cc1f1aa37ff4bcbb109; no amendment).
- **Code:** `scnfet/cg.py` (solver; `scnfet/tests/test_cg.py`) and `scripts/t1_cg.py` (modes `table`, `cache`, `tob`).
  Figure and summary: `scripts/t1_figs.py` → `figures/Fig7_Ion_W.png`, `summary.md`.
- **Outputs:** `cg_table.json` (C_G for every width and DFT point, main geometries and brackets); `single/` and
  `gaa/` (one JSON per width, plus `dft.json`); `summary.md` (every width, both EOTs); logs `cache.log`,
  `tob_single.log`, `tob_gaa.log`. Band cache: `fits/t1/cache/*.npz`.
- **Resources:**
  - band cache: one process, peak 2.04 GB, 19 min;
  - ToB: single gate and GAA in parallel, 0.29 GB each, 1 h 58 min.
- **Cache check:** the cached bands reproduce the parallel-plate I_ON stored in `fits/t1/rows` exactly (relative
  difference < 1e-9) at the first width of each family.
- **Label:** "ballistic, ideal-electrostatics upper bound; geometric C_G (Laplace, HfO2 ε 25)". α_G = 1 and α_D = 0
  stay, so SS is 59.5 mV/dec everywhere, as before. Only C_Σ = C_G changes.

**Decided (16 Nov), after the bare-ribbon numbers:**
- **Parallel plate carries the paper's design rule**, as the convention of the HfS2/QUDEN benchmark.
- **This Laplace result is reported beside it as the other bound,** with the IRDS verdict both ways.
- **Figure 7 shades the range between the two.** Parallel plate is the dense-array limit of the single-gate stack;
  the HfO2 isolated ribbon gives the most fringing.
  - The GAA isolated ribbon is a dotted line above the range. Its own dense-array limit would be 2 × C_pp, which is
    not computed.
  - The registered brackets stay C_G-only.

## Capacitance (C_G / C_pp)
| | EOT 1 nm (HfO2 6.41 nm) | EOT 0.6 nm (HfO2 3.85 nm) |
|---|---|---|
| single gate, HfO2 encapsulated (main) | 3.1 – 10.6 | 2.4 – 7.4 |
| GAA, HfO2 (main) | 4.0 – 12.6 | 3.3 – 9.0 |
| bracket: single gate, vacuum above | 1.9 – 5.8 | 1.6 – 4.1 |
| bracket: single gate, SiO2-equivalent (t = EOT) | 1.5 – 3.2 | 1.3 – 2.5 |
| bracket: GAA, SiO2-equivalent | 2.4 – 4.4 | 2.3 – 3.6 |

The largest values are at 1 nm, where the oxide is thicker than the ribbon is wide.

## I_ON at EOT 0.6 nm (mA/μm; HSE06 building blocks)
Each cell lists three values: parallel plate / single gate / GAA. Every width and both EOTs are in `summary.md`.

| | [100] bare (1.09–5.93 nm) | [110] bare (1.07–6.04 nm) | [110] H, in-sample (1.39–6.36 nm) |
|---|---|---|---|
| n-FET | 2.76–3.80 / 5.94–16.84 / 6.73–18.36 | 1.34–1.61 / 3.98–10.66 / 5.03–13.15 | 2.78–3.82 / 5.77–13.58 / 6.57–14.85 |
| p-FET | 3.07–3.67 / 7.13–18.10 / 9.05–19.52 | 1.53–2.41 / 5.11–11.56 / 6.44–14.16 | 1.43–1.82 / 4.28–9.86 / 5.31–12.31 |
| IRDS, parallel plate | n and p at every width | n 0/17; p 15/17 | n 17/17; p 0/17 |
| IRDS, single and GAA | n and p at every width | n and p at every width | n and p at every width |

**Calibration DFT-Wannier ribbons (relaxed W):**

| ribbon (variant) | W (nm) | n-FET: pp / single / GAA | p-FET: pp / single / GAA |
|---|---|---|---|
| `r100_3` (HW3) | 1.09 | 2.75 / 16.76 / 18.29 | 3.34 / 18.08 / 19.50 |
| `r110_4` (HW1) | 1.07 | 1.31 / 10.49 / 12.98 | 1.48 / 10.92 / 13.32 |
| `r100_3h` (HW2) | 1.32 | 2.87 / 14.48 / 15.79 | 0.81 / 5.16 / 6.42 |
| `r110_4h` (HW2) | 1.39 | 2.66 / 13.35 / 14.65 | 1.42 / 9.84 / 12.29 |

The building blocks reproduce these at the calibration widths under the Laplace C_G as well (for example, [100]
n = 3 single gate: 16.84 / 18.10 against 16.76 / 18.08).

## What it shows
1. **The IRDS verdict flips everywhere.** With the Laplace C_G, every family and every DFT point clears IRDS for n
   and p. The lowest value is [110] bare n at 6 nm, single gate: 3.98 mA/μm, 2.0 × the target. That includes the two
   failures of the parallel-plate design rule, [110] bare n and [110] H p, and the [100] H p-FET (0.81 → 5.16).
2. **The ranking survives, but the edge-band penalties shrink.** The table gives I_ON relative to [100] at equal W
   (1.45–5.9 nm, EOT 0.6):

   | | parallel plate | single gate | GAA |
   |---|---|---|---|
   | [110] bare, n-FET | 0.37–0.47 | 0.63–0.75 | 0.72–0.87 |
   | [110] bare, p-FET | 0.51–0.78 | 0.67–0.72 | 0.66–0.77 |
   | [110] H, n-FET | 0.93–1.01 | 0.96–1.00 | 0.97–1.00 |
   | [110] H, p-FET | 0.40–0.58 | 0.62–0.67 | 0.61–0.75 |

   [100] is the best orientation for both carriers under every C_G. The [110] edge bands always cost current, but
   with a geometric high-k gate they no longer cost IRDS compliance.
3. **The regime changes.** Here is how I_ON at EOT 1 compares with I_ON at EOT 0.6, for the n-FET:
   - **Parallel plate:** 0.56–0.69, close to the 0.6 of a capacitance-limited device.
   - **Laplace, [100] and [110] H:** 0.89–0.94. These ribbons approach the quantum-capacitance limit.
   - **Laplace, [110] bare:** 0.80–0.85. The heavy Sc-s edge band (m ≈ 0.9–1.0) stays further from that limit.

   **My reading (interpretation, not a registered test):** this is why the [110] penalty shrinks.
   - The edge band gains more from the extra charge. With the single gate, I_ON rises 2.5–7.9×, which is 1.01–1.12 ×
     C/C_pp (GAA: 0.94–1.12 ×). Its injection velocity rises 1.4–1.5× as E_F climbs the band (n = 4 and 12).
   - [100] is held back by its own quantum capacitance: I_ON rises 1.7–6.1× for C/C_pp of 2.5–7.2 (single gate).
4. **I_ON now falls monotonically with W,** because fringing matters most for the narrowest ribbons. The
   parallel-plate [100] hump (n at 3.3 nm, p at 2.0 nm) is gone.
5. **This is the most-fringing end.**
   - The geometry is an isolated ribbon, with no neighbouring ribbons screening the field.
   - There is no interfacial layer, and t = EOT · 25/3.9 is the thickest physical oxide for the EOT.
   - Real stacks and arrays lie between this and parallel plate. The SiO2-equivalent bracket (1.3–2.5× at EOT 0.6)
     shows how far down a thin-oxide stack moves.
   - Ballistic transport and ideal electrostatics stay as they were. The Fröhlich B-correction and the WKB guard
     are still to come.
