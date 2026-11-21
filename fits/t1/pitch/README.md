# T1b: pitch-normalised I_ON (laptop, 2026-11-17)

Settings: `fits/PREREG_t1b.md` §1, frozen at 17:16 on 17 Nov before any number (md5 655bbf6f…; no amendment).
`scripts/t1_pitch.py` was hashed at 17:17 before the first run, and re-hashed at 17:18 after a table-format fix
(`pitch.json` byte-identical) in `fits/PREREG_t1b.code.md5`.
- **What it is:** post-processing of `fits/t1/corr/`. I_fp = I_ON × W / (W + s), at EOT 0.6, × B (2D Fröhlich in
  vacuum).
- **Spacings (decision, 17 Nov):**
  - s = 2.0 nm, the gap in our DFT cells, the densest array the Hamiltonians cover. GAA doesn't fit there.
  - s = 7.69 nm = 2 t_HfO2, where the oxide shells touch (the smallest GAA pitch).
- **Fill factor W/P:** 0.35 → 0.75 at s = 2 nm, and 0.12 → 0.44 at s = 7.7 nm, over 1.1–6 nm.
- **Outputs:**
  - `pitch.json`, `summary.md` (every width, plus the DFT calibration ribbons);
  - `figures/Fig8_Ion_pitch.png`, with three columns: per W, then per footprint at s = 2 nm and at s = 7.7 nm.

## IRDS pass counts, × B (ribbons passing / total)
| | per W: pp / single / GAA | s = 2 nm: pp / single | s = 7.7 nm: pp / single / GAA |
|---|---|---|---|
| [100] bare n | 10 / 12 / 12 | **0 / 12** | 0 / 0 / 0 |
| [100] bare p | 9 / 12 / 12 | **0 / 12** | 0 / 1 / **8** |
| [110] bare n | 0 / 4 / 8 | **0 / 0** | 0 / 0 / 0 |
| [110] bare p | 0 / 16 / 16 | **0 / 16** | 0 / 0 / 0 |
| [110] H n | 15 / 17 / 17 | **0 / 17** | 0 / 0 / 0 |
| [110] H p | 0 / 6 / 11 | **0 / 0** | 0 / 0 / 0 |

Ballistic (before B), the same columns: `summary.md`.

## Per-footprint I_ON × B (mA/μm), lowest to highest width
| | s = 2: pp | s = 2: single | s = 7.7: pp | s = 7.7: single | s = 7.7: GAA |
|---|---|---|---|---|---|
| [100] n | 0.54–1.60 | 2.72–3.27 | 0.19–0.93 | 1.15–1.63 | 1.26–1.85 |
| [100] p | 0.74–1.47 | 3.41–4.02 | 0.26–0.85 | 1.41–1.99 | 1.52–2.52 |
| [110] bare n | 0.14–0.35 | 0.87–1.14 | 0.05–0.21 | 0.40–0.51 | 0.49–0.65 |
| [110] bare p | 0.35–1.07 | 2.03–2.32 | 0.13–0.62 | 0.76–1.33 | 0.92–1.67 |
| [110] H n | 0.61–1.66 | 2.70–3.13 | 0.23–0.99 | 1.11–1.66 | 1.22–1.89 |
| [110] H p | 0.18–0.44 | 1.03–1.27 | 0.07–0.26 | 0.47–0.61 | 0.59–0.76 |

([110] bare p leaves out n = 4, the flat valence band, which gives 0.02 / 0.007 under parallel plate.)

## What it shows
1. **Per footprint, the family ranking is the same as per W.** The fill factor is common to all families at a
   given W, so the edge-band carriers stay the weakest: [110] bare n and [110] H p, at about a quarter to a third of [100].
2. **The IRDS verdict per footprint depends on the array electrostatics, and the bracket spans everything.**
   - **Parallel plate × B:** no ribbon passes at either spacing.
   - **Isolated single gate at s = 2 nm:** the design rule holds cleanly at every width. [100] n and p, [110] bare p
     and [110] H n pass (all widths, except the [110] bare flat-band point). The two edge-band carriers fail
     everywhere.
   - **At the GAA spacing:** only the [100] p-FET passes, with GAA from 2.85 nm (8 of 12) and the single gate at
     5.9 nm.
3. **The array C_G was computed afterwards (17 Nov, `fits/t1/cg_array/README.md`; approved).** It closes the
   bracket.
   - Per W, the design rule holds at both spacings.
   - Per footprint at s = 2 nm, the light-band carriers sit around the IRDS line (1.5–2.15 mA/μm, 0.77–1.09 of it) and the edge-band
     carriers 3–4× below. At 7.7 nm nothing passes with a single gate.
   - *(Originally written here: "The true array C_G at these spacings is not computed … the isolated-ribbon end is
     optimistic there." The array value confirmed that expectation.)*
4. **Narrow ribbons pay their fill factor.** Per footprint under parallel plate, I_ON rises with W for every family,
   while per W it falls. The sheet limit (W → ∞) has a fill factor of 1.
5. **[110] bare p at s = 2 nm, single gate** sits just above the line (2.03–2.32). Its pass is fragile.

## Limits to state
- **Every limit of `fits/t1/corr/README.md` carries over:** ballistic ToB, ideal electrostatics, α > 1 on the heavy
  edge bands, vacuum ε12, remote HfO2 phonons left out.
- **No gate metal between ribbons at s_GAA**, so it is a lower bound on the real GAA spacing.
- **The headline verdict stays per W, as registered.** Per footprint is the referee-requested comparison with a
  sheet of the same footprint.
