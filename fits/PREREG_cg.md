# Geometric gate capacitance (Laplace C_G), fixed before any device number with it

Written 2026-11-16, before any I_ON was computed with a Laplace C_G. Go-ahead: 16 Nov.
- **What it is:** plan D6 / referee item K3, "geometric C_G(W) by 2D Laplace, single-gate and GAA, beside
  parallel-plate".
- **Status:** a labelled correction to the T1 device numbers. Parallel-plate stays the signed-off convention, and
  the HfS2 benchmark never uses this.
- **Freezing:** `fits/PREREG_cg.md5`; any later change is a dated amendment below.

## Geometry (solver: `scnfet/cg.py`)
- **Ribbon:** a perfect conductor of width W (edge to edge, the same W as in T1, from the copy-step geometry) and
  thickness d = 0.217 nm (the two planes, nuclei to nuclei). All families use the same d. H atoms count only through
  W, as in the convention.
- **Gate dielectric: HfO2, ε_r = 25, with physical thickness t = EOT × 25 / 3.9.** That gives 6.41 nm at EOT 1 nm and
  3.85 nm at EOT 0.6 nm.
  - The EOT stays referred to SiO2 (3.9), as the signed-off convention defines it, so the wide-ribbon limit equals the
    parallel-plate value.
  - Why HfO2 rather than SiO2 of thickness EOT: an EOT of 0.6–1 nm is only reached with high-k stacks (IRDS). The
    physical thickness, not the EOT, sets the fringing.
- **Single gate (main):** a planar gate. The ribbon's bottom face sits at height t, and the ribbon is encapsulated:
  HfO2 everywhere above the gate.
- **GAA:** a rectangular gate at distance t from every face of the ribbon, HfO2 inside.
- **Descriptive brackets, C_G only, with no device numbers:**
  - single gate with vacuum above the ribbon (HfO2 only between gate and ribbon);
  - SiO2-equivalent (ε 3.9, t = EOT) for both geometries, the least-fringing bound.

## Numerics (in-sample checks before freezing)
- **Grid:** finite volumes; step h = 0.02 nm at every corner and interface, growing ×1.15 towards segment middles
  and outwards. Single gate: zero-flux outer boundary at 40 × max(W, t).
- **Microstrip check:** a zero-thickness strip in a homogeneous medium matches the Hammerstad–Jensen closed form
  within 0.1–1.4 % (W/h 0.5–6), about that formula's own accuracy.
- **Convergence:** h 0.02 → 0.01 changes C by 0.05 %; extent 40 → 80 × max(W, t) changes it by 0.04 %.
- **Wide limits:** the single gate at W/t = 200 gives 1.024 × parallel plate; GAA gives 1.005 × twice parallel plate.

## Device numbers with the Laplace C_G
- **Same as the parallel-plate T1 run:** the ToB code and conventions; the HSE building-block bands, recomputed by
  the same T1 code and cached; the I_OFF pin; and the doping.
  - Only C_Σ = C_G changes, per unit length, from the geometry above at each family's W.
  - The I_ON per W keeps the same W.
- **Families:** [100] bare, [110] bare and [110] H at every T1 width, plus the four calibration DFT device points
  (at their relaxed W).
- **Geometries:** single gate and GAA, n and p, EOT 1 and 0.6 nm.
- **Reported:** I_ON, SS and the IRDS verdict (≥ 1.979 mA/μm at EOT 0.6), beside the parallel-plate results. This is
  the plan's "IRDS verdict reported both ways".
- **Label:** "ballistic, ideal-electrostatics upper bound; geometric C_G (Laplace, HfO2 ε 25)". Ideal electrostatics
  (α_G 1, α_D 0) stays; only the capacitance changes.
- **Report everything; nothing is tuned.**
