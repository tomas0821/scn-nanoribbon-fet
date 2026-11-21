# Pre-registered: the SK hold-out tests, claims gates (b) and (d)

Written 2026-11-16, **before any SK prediction for the slab or the hold-out ribbons was computed.** The DFT and
Wannier references were already known (A3, A5 READMEs); the SK side was not. Go-ahead: 16 Nov.
The recipe choices stay within the wording of plan §3 D3, D4 and D7; they are frozen by the md5 in
`fits/PREREG_sk_holdouts.md5`, which `scripts/gate_bcd_sk_holdouts.py run` checks before it computes anything. Any
later change is a dated amendment below, with both results.

Already computed before freezing (no hold-out data): the self-test and gate (c), both in-sample (see the end).

## Models (fixed; no refit)
- **PBE twin:** `fits/A2_joint/sk_joint_PBE.json`. **SK-HSE:** `fits/A2_joint/sk_joint_HSE06.json`. These are the
  pre-registered joint bulk + sheet fits as they stand. Gate (a) failed on them (2026-11-11; LR1 → DFT-Wannier device
  Hamiltonians); SK remains the model that interpolates across widths, and these tests measure its transfer errors.
- **Form:** as `fits/PREREG_sk_bulk_sheet.md`. Shells `rocksalt_shells(4.501)` with `DEFAULT_ETA`, the three-centre
  term (η 3.5), on-site Ep / Et2g / Eeg with t2g/e_g on the cubic axes (`Structure.frame`).
- **Geometry:** the DFT-relaxed structures the references were computed at: the slab from A3-1
  (`make_a3.relaxed()`), the hold-outs from A5-1 (`make_a5.relaxed`), the calibration ribbons from A4a
  (`make_a4b.geometry`). Hoppings are SK with Harrison scaling at those distances. The T1 step that copies the
  calibration edge relaxation onto wider ribbons (D4) is not tested here; the A5-1 edges relax like the A4a ones
  (edge-atom environments agree to 0.01 Å, self-test 3).
- **Surface on-site blocks:** the sheet's fitted C4v shifts on every atom of a surface plane: all atoms of the
  two-plane ribbons, and the two outer planes of the slab. The slab's middle plane has none (bulk on-site). The shift
  patterns are rotated from the sheet's cell frame into each structure's frame (self-test 1).
- **Edge on-site blocks** (D4, "orbital-resolved edge blocks read from W90 H(R)"):
  - For each atom of the outermost host row of each edge:
    B = H_W(e, R = 0) − ⟨H_W(c, R = 0)⟩_c, added to that atom's SK on-site block.
  - H_W is the calibration ribbon's Wannier Hamiltonian with the same direction and edge chemistry: [100] bare
    n = 3 for [100] n = 6, and [110] bare n = 4 for [110] n = 8.
  - e is the calibration edge atom with the same side, species and plane.
  - c runs over the atoms of the same species in the calibration ribbon's two centre rows.
  - The blocks are full d (5 × 5) or p (3 × 3) matrices in the model's orbital order. Wannier90 orbitals are mapped by
    `wf_map.json`, with signs +1 (the A1 check).
  - Nothing else depends on the edge.
  - **Sources:** the variants the A4b decisions named for the SK blocks. PBE twin: A4b-1 W2 for both ribbons.
    SK-HSE: A4b-2 `r100_3` HW2 and `r110_4` HW1 (amendment 2).
  - **Orbitals outside the SK basis are dropped:** the [110] edge Sc s WFs.
    - So the SK model has no Sc-s edge band. In bare [110] the lowest PBE-Wannier conduction band is that band
      (LR3), and the registered gap, valley and CB-mass tests compare it with SK's lowest conduction band.
    - A descriptive line also gives the Wannier values without the edge-s bands (bands with edge-s weight ≥ 0.5
      left out).
- **No energy reference** is needed: every metric is a difference.

## Metrics (the same code for SK and Wannier)
**Slab (gate b), Γ–X̄–M̄–Γ with 200 points per segment (`parse_a3.frac_path`):**
- gap: CB minimum − VB maximum on the path;
- E_c(M̄) − E_c(Γ) and E_c(X̄) − E_c(Γ);
- VB and CB masses along Γ–X̄ and Γ–M̄ by `parse_a4b.edge_mass` (parabolic, the points within 30 meV of the extremum,
  at least 3), applied to each model on its own; VB−1 is descriptive.
- **Reference:** the PBE-Wannier slab W2 (A3 decision).
- **HSE:** the SK-HSE path gap against the HSE grid gap, 1.1055 eV (direct at Γ).
- **Descriptive:** E_c(X̄) − E_c(Γ) also against DFT, because W2 misses it by 0.49 eV (the claims-gate flag).

**Ribbons (gate d), on a 601-point path Γ → X (k_frac 0 … ½):**
- gap: lowest CB minimum − highest VB maximum, with n_val = 3·N_N valence bands (= nocc − nexcl, asserted);
- **valley offset** E_c(X) − E_c(Γ) of the lowest CB at the zone boundary and the zone centre. In [110] this is the
  sheet's M̄ valley, which folds to X. In [100] the M̄ valley folds onto Γ and X carries the X̄ valley, so there this
  tests the X̄-derived valley;
- VB and CB band-edge masses along the ribbon, by `edge_mass` as above;
- **I_ON** from `scnfet.tob` with the signed-off conventions:
  - n and p, EOT 1 and 0.6 nm, VDD 0.7 V, I_OFF 10 nA/μm, 300 K;
  - parallel-plate C_G, with W = `structures.width` of the relaxed hold-out (nuclei, edge to edge) in C_G and I/W;
  - doping 0.01 per atom;
  - bands on a uniform 2000-k grid, with a 4000-k check on the Wannier side.
- **IRDS verdict:** at EOT 0.6 nm, I_ON ≥ 1.979 mA/μm, for n and p separately. It passes if SK and Wannier give the
  same yes/no for both. This is the I_ON comparison only. The PBE twin gaps are below the 0.9 eV BTBT criterion, which
  is a T1 matter, not part of a twin-against-twin verdict.
- **Reference:** the PBE-Wannier hold-outs W2 (A5-2 decision, H4 ≤ 3.6 meV).
- **HSE:** the SK-HSE path gap against the A5-3 HSE grid gaps: 0.921 eV (`r110_8`) and 1.629 eV (`r100_6`), both
  direct at Γ.

## Pass criteria (as registered: plan D7; A3 README "Hold-out criteria")
- **PBE twin vs PBE-Wannier, ribbons:**
  - gap ±50 meV;
  - valley offset ±50 meV;
  - masses ±10 % (VB and CB);
  - I_ON ±15 % (all four carrier/EOT values);
  - the same IRDS verdict.
- **PBE twin vs PBE-Wannier, slab:**
  - gap ±50 meV;
  - both valley offsets ±50 meV;
  - the four band-edge masses ±10 %.
- **SK-HSE gaps:** ±80 meV on each hold-out.
- A hold-out passes only if every item passes.
- **Failure → plan LR2:** no refitting; reported as a transferability limit. The device-Hamiltonian route that
  follows is decided at the claims gate (the plan's LR2: building-block Wannier ribbons for widths ≤ 6 nm).

## Descriptive (no pass/fail)
- The same metrics for SK **without** edge blocks, to show what the blocks do.
- **Edge-block transfer (the D4 drift check):** the hold-out's own W2 edge blocks, by the same formula, against the
  calibration blocks (largest element difference). This uses hold-out data and enters no SK prediction.

## Already computed before freezing (in-sample; `fits/holdouts/selftest.json`, `gate_c.json`)
- **Self-test:**
  1. A [100]-aligned sheet supercell with the rotated shifts reproduces the primitive sheet at the folded k to
     1e-14 eV. The unrotated shifts miss by 0.09 eV (PBE) and 1.55 eV (HSE), so the test detects a wrong frame.
  2. The constant edge blocks land on the right matrix block.
  3. The hold-out edge atoms map to calibration edge atoms with the same environment (≤ 0.010 Å).
  4. The slab has 4 surface-plane atoms of 6.
- **Gate (c), SK + its own edge blocks vs the Wannier Hamiltonian they came from:**

| Ribbon | Functional | Wannier gap | SK + blocks | SK, no blocks |
|---|---|---|---|---|
| [100] n = 3 | PBE | 0.915 eV | + 42 meV | + 86 meV |
| [100] n = 3 | HSE | 1.964 eV | + 82 meV | + 79 meV |
| [110] n = 4 | PBE | 0.419 eV | + 1,040 meV | + 1,231 meV |
| [110] n = 4 | HSE | 1.053 eV | + 1,314 meV | + 1,664 meV |

  - **[100] n = 3:** CB mass within 1.1 %; VB mass − 16 % (PBE) and − 32 % (HSE), carried over from gate (a).
  - **[110] n = 4:** the gap is set by the Sc-s edge CB, which SK lacks. Without the edge-s bands the Wannier gaps are
    1.206 eV (PBE) and 2.100 eV (HSE), so SK + blocks is still + 253 / + 267 meV high on the interior conduction band.
  - **So the [110] n = 8 test is expected to fail on the gap.** It is registered as it stands: no edge s is added to
    the SK form after seeing this. That would be a new model, and a separate decision.
