# Pre-registered: building-block Wannier ribbons on the PBE hold-outs (plan LR2)

Written 2026-11-16 (16:20), **before any building-block Hamiltonian of a hold-out was assembled.** The hold-out
Wannier references were already known (A5 README; gate (d), `fits/holdouts/README.md`); the building-block side was
not. Go-ahead: the claims-gate decision (16 Nov). The recipe choices stay within plan
LR2 and D7; they rest only on the sheet and the calibration ribbons (in-sample checks below). They are frozen by the
md5 in `fits/PREREG_bb_holdouts.md5`, which `scripts/gate_bb_holdouts.py run` checks before it reads any hold-out
Hamiltonian. Any later change is a dated amendment below, with both results.

**This is the second test on the same two hold-outs.** The first was the SK test (gate (d), failed, 16 Nov). The
hold-outs are never refitted, and nothing here is fitted at all.

**Known before freezing that bears on the recipe:** gate (d)'s descriptive drift check, which used hold-out data:
the hold-outs' edge on-site blocks, referred to their own centre rows, differ from the calibration ribbons' by
≤ 48 meV (the [110] Sc d nearly uniformly by ≈ + 40 meV). That was known when the energy alignment below was chosen.

## Recipe (code: `scnfet/bb.py`, driver `scripts/gate_bb_holdouts.py`)
- **Method:** building blocks after Lee, Nardelli & Marzari (PRL 95, 076804), plan D7 and LR2: the calibration
  ribbons' edge regions joined to the sheet's interior, laptop only, existing DFT only.
- **Interior source:** the A2 sheet Wannier Hamiltonians `dft/A2/w90/{pbe_k16,hse_k16}` (16 projection-only WFs,
  16 × 16 k, A2's window rule).
- **Edge source:** the calibration ribbon with the same direction and edge chemistry, `r100_3` for `r100_6` and
  `r110_4` for `r110_8`, in **variant W1 (PBE) / HW1 (HSE)**: A2's window rule as run (A4b README, W1), the same rule
  the sheet's Wannier functions were made with. W1 passed P1–P5 for both ribbons; HW1 passed H2.
  - Why not W2, the SK-block variant: the in-sample check below shows that W2's Sc z² Wannier functions differ from
    the sheet's (the sheet's sit 0.94 Å out into the vacuum, W2's 0.1–0.3 Å), so W2 and sheet elements disagree by
    0.2–0.9 eV on z². With W1, the [100] centre-row elements agree with the sheet's within 40 meV.
  - The hold-out reference stays W2 (A5-2 decision). The bands are compared, not the matrix elements, and the
    hold-outs' W1 and W2 agree within H4 ≤ 3.6 meV over the device window.
- **Basis:** Sc d + N p on every atom, and Sc s on the bare [110] edge Sc (taken with the calibration ribbon's edge
  s), in the Wannier90 real-harmonic order. That is the reference's basis.
- **Assembly:**
  - Atoms are matched on the ideal cut by row counted from each edge, plane, and position along x (mod the period).
  - An element is copied from the calibration ribbon when both atoms lie within the calibration ribbon's half
    width of the same edge: N_C / 2 rows, which is 3 rows for [100] and 4 for [110].
  - Every other element is the sheet's element for the same species, planes and in-plane bond vector, with the
    orbitals rotated from the sheet's axes into the ribbon's: Q = R45 F_Tᵀ, the identity for [110].
  - Bond vectors come from the ideal cut. The calibration edge relaxation enters through its matrix elements: this
    is the D4 step that copies the calibration edge relaxation onto wider ribbons.
  - Elements with |d| > 20 Å are zero. In-sample, 12 Å moved the `r110_4` gap by 21 meV, and 20 Å moves it by
    ≤ 1.2 meV.
- **Energy alignment:**
  - Each source is first referred to its own vacuum level: the A2 / A4b E_vac values.
  - The sheet is then shifted by one constant: the mean on-site diagonal of the calibration ribbon's two centre rows
    (all atoms and orbitals) minus the sheet's. This is Lee/Nardelli/Marzari's rule of matching on-site energies
    where both sources describe the same environment.
  - In-sample, the narrow [100] ribbon's centre sits a uniform ≈ 0.21 eV below the vacuum-aligned sheet, with Sc and
    N alike. Vacuum alignment alone would put a 0.2 eV step between the edge region and the interior.
- **Nothing is fitted, and nothing depends on the hold-outs** apart from their ideal cut.

## Metrics and pass criteria (plan D7, exactly as gate (d): `fits/PREREG_sk_holdouts.md`)
- **Shared with gate (d):** the same code (`gate_bcd_sk_holdouts.ribbon_metrics`, `compare`, `ion`), on the
  601-point Γ → X path and the 2000-k ToB grid, with the gate sweep ±4 V / 381 points.
- **Reference:** the PBE-Wannier hold-out W2.
- **BB PBE against W2:**
  - gap ±50 meV;
  - valley offset E_c(X) − E_c(Γ) ±50 meV;
  - VB and CB masses ±10 %;
  - I_ON ±15 % for n and p at EOT 1 and 0.6 nm;
  - the same IRDS verdict (≥ 1.979 mA/μm at EOT 0.6 nm, n and p).
- **BB HSE06:** the gap against the A5-3 HSE grid gap, 0.921 eV (`r110_8`) and 1.629 eV (`r100_6`), ±80 meV.
- **A hold-out passes only if every item passes.**
- **Descriptive (no pass/fail):**
  - BB PBE with vacuum alignment only (sheet shift 0);
  - BB PBE with W2 edges (the SK-block variant, gauge-mismatched);
  - for [110], the gap and valley without the edge-s bands.

## What follows (decision of 16 Nov, read per direction)
- **A direction whose bare hold-out passes:** building blocks become that direction's T1 device Hamiltonians at
  1–6 nm. Its H-terminated ribbon gets a building-block curve labelled "validated on bare edges only".
- **A direction whose hold-out fails:** the decided fallback.
  - [100]: SK with stated errors, plus the DFT-Wannier points.
  - [110]: DFT-Wannier points only, with no width curve, and its H ribbon likewise.
- Reading the decision per direction, rather than all-or-nothing, is this document's interpretation; it can be
  overruled.
- No recipe change after the run except as a dated amendment, with both results.

## Already computed before freezing (in-sample; `fits/bb/selftest.json`, `fits/bb/calib.json`)
- **Self-test:**
  1. A [100]-aligned conventional sheet cell built from rotated sheet elements reproduces the primitive sheet at the
     folded k to 1e-14 eV (PBE and HSE).
  2. `orb_rot` is orthogonal and composes, and the 90° signs match the real harmonics (`scnfet/tests/test_bb.py`;
     58 tests pass).
  3. Row counts: `r100_6` 12 rows / `r100_3` 6; `r110_8` 16 / `r110_4` 8.
- **Calibration ribbon rebuilt from its own halves plus the sheet** (h = N_C / 2: only cross-half pairs from the
  sheet), against its own Wannier bands:

| Ribbon, variant | sheet shift | gap | valley | m_VB | m_CB | device-window dev | centre-row hopping vs sheet (max) |
|---|---|---|---|---|---|---|---|
| `r100_3` PBE W1 (recipe) | − 0.212 eV | − 5 meV | 0 | − 0.1 % | + 0.8 % | 13 meV | 39 meV |
| `r100_3` PBE W2 | − 0.060 | − 5 | + 4 | + 0.2 % | − 6.0 % | 50 | 334 (z²) |
| `r100_3` HSE HW1 (recipe) | − 0.243 | − 12 | + 2 | − 0.6 % | − 0.3 % | 40 | 100 |
| `r100_3` HSE HW2 | − 0.109 | − 9 | + 8 | − 0.2 % | − 3.5 % | 34 | 266 (z²) |
| `r110_4` PBE W1 (recipe) | + 0.139 | − 27 | − 18 | **− 39.7 %** | − 4.7 % | 46 | 149 |
| `r110_4` PBE W2 | + 0.257 | − 27 | − 21 | − 33.1 % | − 6.1 % | 61 | 320 (z²) |
| `r110_4` HSE HW1 (recipe) | + 0.196 | − 24 | − 2 | **+ 87.9 %** | − 6.2 % | 54 | 198 |

- **[100]:** the sheet and the W1 ribbon join cleanly. The rebuild reproduces the ribbon within 5–12 meV in the gap
  and 1 % in the masses.
- **[110]: the narrow ribbon's centre is not sheet-like even in the shared gauge.** Its centre rows sit + 0.11 to
  + 0.20 eV off the sheet, unevenly by species and row, and its hoppings are up to 149 meV off. Putting sheet
  elements only between its two halves already moves the VB mass by − 40 % (PBE) and + 88 % (HSE).
- **So the [110] n = 8 test is expected to fail on the hole mass,** unless its interior is more sheet-like than
  n = 4's. It is registered as it stands.
- Changing the edge depth or alignment after this would be a new recipe, and a separate decision.
