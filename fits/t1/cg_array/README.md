# Gate capacitance of a ribbon array (laptop, 2026-11-17)

Settings: `fits/PREREG_array.md`, frozen at 18:12 on 17 Nov (md5 d27cea0c…; no amendment). Code md5s are in
`fits/PREREG_array.code.md5`, recorded before the run. `scripts/t1_array_post.py` was hashed at 18:27 while the ToB
ran.
- **Solver:** `scnfet/cg_array.py`, one period with mirror planes, using `cg.py`'s solver. Tests:
  `scnfet/tests/test_cg_array.py` (71 tests pass in total).
- **Run:** `fits/t1/cg_array/run_tob.sh`, two processes at EOT 0.6 from the T1 cache. 1 h 03 min each, peak 0.26 GB.
- **Outputs:**
  - `cg_array_table.json` (C_G);
  - `s2/`, `sGAA/` (ToB per ribbon, plus `dft.json`);
  - `array.json` and `summary.md` (× B, per W and per footprint).
- **Check:** for all 92 ribbon-carrier pairs, I_ON(pp) ≤ I_ON(array, 2 nm) ≤ I_ON(array, 7.7 nm) ≤ I_ON(isolated).

## Capacitance
| | s = 2 nm | s = 7.69 nm |
|---|---|---|
| C / C_pp | 2.7 → 1.3 (1.1 → 6 nm) | 5.4 → 1.8 |
| C / C_isolated | 0.37 → 0.54 | 0.74 → 0.76 |

At s = 2 nm, below the 3.85 nm oxide, neighbours screen most of the fringing field. GAA needs no array solve: at
s ≥ 2t each gate box encloses its ribbon.

## IRDS × B, EOT 0.6 (ribbons passing / total)
| | per W, s = 2 | per W, s = 7.7 | per footprint, s = 2 | per footprint, s = 7.7 |
|---|---|---|---|---|
| [100] n | 12 / 12 | 12 / 12 | 2 / 12 (2.85, 3.29 nm) | 0 |
| [100] p | 12 / 12 | 12 / 12 | 5 / 12 (≤ 2.85 nm) | 0 |
| [110] bare n | **0** / 17 | **1** / 17 (1.07 nm) | 0 | 0 |
| [110] bare p | 4 / 17 (1.7–2.6 nm) | 16 / 17 | 0 | 0 |
| [110] H n | 17 / 17 | 17 / 17 | 5 / 17 | 0 |
| [110] H p | **0** / 17 | **1** / 17 (1.39 nm) | 0 | 0 |

Per-footprint I_ON × B at s = 2 nm (mA/μm):
- light-band carriers: [100] n 1.58–1.99, [100] p 1.91–2.15, [110] H n 1.52–2.02;
- [110] bare p: 0.79–1.37;
- edge-band carriers: [110] bare n 0.39–0.47, [110] H p 0.44–0.58.

## What it shows
1. **With the array's own capacitance, the design rule holds per W at both spacings.**
   - [100] passes for both carriers at every width, and [110] H n at every width.
   - The two edge-band carriers fail at every width at s = 2 nm, and at all but the narrowest at s = 7.7 nm.
   - The bracket of Figure 8 is closed: the array sits between parallel plate and the isolated ribbon, as it must.
2. **Per footprint, a dense array puts the light-band carriers around the IRDS line** after backscattering, at
   0.77–1.09 of it, while the edge-band carriers sit a factor 3–4 below it. At the GAA spacing nothing passes per footprint
   with a single gate. With GAA, only the [100] p-FET passes, from 2.85 nm (`fits/t1/pitch/`).
3. **The headline verdict stays per W, as registered.** The array result makes it stronger. Under the array
   capacitance the rule no longer depends on choosing between parallel plate and the isolated ribbon.
