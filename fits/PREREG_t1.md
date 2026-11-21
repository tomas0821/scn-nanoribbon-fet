# T1 settings, fixed before the width sweep: building-block device Hamiltonians, 1–6 nm

Written 2026-11-16, before any T1 width result. Go-ahead: 16 Nov, after the building-block ribbons passed
both hold-outs (`fits/bb/README.md`). This is production rather than a test: it fixes what is computed and how, so
nothing is chosen after seeing the curves. It is frozen by `fits/PREREG_t1.md5`, and any later change is a dated
amendment below.

## Device Hamiltonians
- **Production: HSE06 building blocks**, by the frozen recipe of `fits/PREREG_bb_holdouts.md`:
  - the A2 HSE sheet interior;
  - the calibration ribbon's HW1 edges: `r100_3`, `r110_4`, `r100_3h` and `r110_4h`;
  - the centre-row alignment and the 20 Å cutoff.
- **The PBE building blocks** (W1 edges) are computed for bands only, to sit beside the PBE-Wannier DFT points.
- **Families and widths:** every integer n in the plan's range.
  - bare [100] and [100] H, n = 3–14 (1.1–5.9 nm);
  - bare [110] and [110] H, n = 4–20 (1.1–6.1 nm).
- **Geometry:** the ideal cut, with the matched calibration atoms' relaxation displacements on the edge-region
  atoms (D4's copy step) and the rigid shift removed. W is measured from this geometry, edge to edge with H included.
  This is checked against the two relaxed hold-out widths.

## H-terminated ribbons (no hold-out tests them)
- **First, an in-sample check:** each H calibration ribbon is rebuilt from its own halves plus the sheet, as in
  `fits/bb/calib.json`, in PBE (W1) and HSE (HW1). The rebuild is compared with that ribbon's own bands.
- **Pass:** gap ± 50 meV, valley ± 50 meV and both masses ± 10 %, in both functionals. The ribbon then gets a
  building-block curve labelled "validated on bare edges only; in-sample for H".
- **Fail:** no curve for that H family. It is reported at its DFT point only, and the failure goes back to the authors.
- This rule is stricter than the bare hold-out evidence would require (the bare [110] in-sample rebuild failed, yet
  its hold-out passed). With no hold-out, the in-sample check is the only evidence for H.

## Metrics (the conventions signed off 11 Nov; code as gates (d) and BB)
- **Bands:** 601-point Γ → X path.
  - gap, and whether it is direct;
  - the valley offset E_c(X) − E_c(Γ) of the lowest CB;
  - the VB and CB band-edge masses (`edge_mass`);
  - for bare [110], the gap and valley without the edge-s bands (bands with edge-s weight ≥ 0.5 left out).
- **ToB:** n and p, EOT 1 and 0.6 nm, VDD 0.7 V, I_OFF 10 nA/μm, 300 K, doping 0.01 per host atom.
  - Parallel-plate C_G with W edge to edge.
  - 2000 k, and a gate sweep of ±4 V in 381 points.
  - Outputs: I_ON, SS, v_inj, and the IRDS verdict (I_ON ≥ 1.979 mA/μm at EOT 0.6 nm).
  - All labelled "ballistic, ideal-electrostatics upper bound; parallel-plate C_G". The Laplace C_G, the Fröhlich
    B-correction and the WKB guard are later T1 items and do not change these numbers.
- **Valley reordering:** a sign change of E_c(X) − E_c(Γ) with W.
  - It is reported as a result only where DFT ribbons bracket it: [100] n = 3–6, [110] n = 4–8.
  - Anywhere else it is a "model prediction".

## Overlays and stated errors
- **DFT points:**
  - HSE grid gaps: A4b-2 for the calibration ribbons, A5-3 for the hold-outs;
  - ToB on the calibration HSE device Hamiltonians as decided in A4b-2: `r100_3` HW3, `r110_4` HW1, `r100_3h` HW2,
    `r110_4h` HW2;
  - PBE-Wannier gaps (W2) beside the PBE building blocks.
- **SK-HSE on bare [100]:** the gate (d) recipe (joint fit, HW2 edge blocks) on the same geometry, as a cross-check.
- **Stated errors:** the building-block validation errors at the hold-outs, BB against W2 and HSE.
  - gap − 8 meV ([100]) / − 32 meV ([110]);
  - HSE gap − 25 / − 10 meV;
  - masses ≤ 0.5 % / ≤ 8 %;
  - I_ON ≤ 1.3 % / ≤ 3.9 %.
- These are the error bars quoted for now. A width-dependent error model is a later item.
