# Pre-registration T1b: pitch-normalised I_ON (Figure 8) and E_g(W) error bands (Figure 4)

Frozen 2026-11-17, before any pitch-normalised number or error band is computed. The settings are the decisions
of 17 Nov. Later changes are dated amendments below; this text is not edited.

Both items are post-processing of numbers that already exist (`fits/t1/rows/`, `fits/t1/corr/`). No new band
structure, ToB or C_G solve is run.

## 1. Pitch-normalised I_ON (plan D6, LK11; Figure 8)
- **Definition:** I_fp = I_ON(W) × W / P, with P = W + s, in mA per μm of array footprint.
  - W is edge to edge, the same W as the per-W normalisation.
  - s is the edge-to-edge gap between neighbouring ribbons.
- **Two spacings (the bracket, as decided):**
  - **s_min = 2.0 nm:** the inter-ribbon gap in our DFT cells (20 Å for [110] bare and H and [100] H, with edge
    bands converged to ≈ 10 meV; bare [100] is converged at 15 Å). It is the densest array the Hamiltonians cover.
  - **s_GAA = 2 t_HfO2 = 2 × EOT × 25 / 3.9 = 7.69 nm at EOT 0.6.** At this spacing the oxide shells of the
    `scnfet/cg.py` stack touch, so it is the smallest GAA pitch. With no gate metal between ribbons, it is a lower
    bound on a real GAA spacing.
- **Inputs:** I_ON at EOT 0.6 from `fits/t1/corr/<fam>_n<n>.json`:
  - ballistic and × B (2D Fröhlich in vacuum, main);
  - parallel plate, single gate and GAA.
- **Gate capacitance at finite pitch:** the per-W bracket is used unchanged.
  - **The bracket is valid at any spacing.** Neighbours at the same potential screen each other's fringing, so an
    array's C_G per ribbon width lies between parallel plate (the s → 0 limit) and the isolated ribbon (s → ∞).
    The ToB I_ON is monotonic in C_G.
  - **GAA is reported at s_GAA only.** At s = 2 nm it doesn't fit geometrically: the oxide shell alone is 3.85 nm
    per side.
  - **The single gate is reported at both spacings.**
- **IRDS:**
  - **The plot:** the 1.979 mA/μm line is drawn on every panel.
  - **Pass counts:** quoted for each spacing, geometry and B setting. Per footprint, an array is compared with a
    sheet device of the same footprint, for which per W and per footprint are the same.
  - **The headline verdict stays per W, as registered.**
- **Outputs:**
  - `scripts/t1_pitch.py` → `fits/t1/pitch/pitch.json` and `summary.md`;
  - `scripts/t1_figs.py pitch` → `figures/Fig8_Ion_pitch.png`.

## 2. E_g(W) error bands (Figure 4)
- **Model (as decided): a constant symmetric band per direction**, equal to the larger of the PBE and HSE
  hold-out gap errors of the building blocks (`fits/bb/README.md`):
  - **[100]:** ± 25 meV (`r100_6`: PBE − 8, HSE − 25 meV);
  - **[110]:** ± 32 meV (`r110_8`: PBE − 32, HSE − 10 meV).
- **Constant in W:** the edge on-site blocks do not drift with width (≤ 0.05 eV, gates (b)–(d)), and the interior
  converges to the sheet, which is DFT by construction. So no growth with W is assumed.
- **The sign is not used:** all four hold-out errors are negative (building-block gaps low). They are reported, not
  used to shift the curves; shifting would be tuning to the hold-outs.
- **The same band goes on both panels** (HSE06 production, PBE twin).
- **[110] H borrows the [110] band** and keeps its "in-sample only" label. [100] H stays a DFT point.
- **Not added to the band, stated in the caption:**
  - the 20 Å vacuum residual, 8–12 meV, on edge-band states only (2026-11-13 decision);
  - the SK-HSE [100] cross-check, whose error against the building blocks is shown directly by its curve.
- **Reference points at W → ∞** (DFT, the thickness axis is not claimed): the two-plane sheet (HSE 1.415 eV,
  PBE 0.404 eV; A2) and the three-plane slab (HSE 1.106 eV, PBE 0.120 eV; A3), as right-margin markers.
- **Plan wording:** plan D4 says "error bands from cross-validation and the fit covariance". That was written for SK
  as the production model. Under LR2 the device Hamiltonians are building blocks with no fit, so the plan gets a
  dated amendment pointing here.

## Amendments
*(none)*
