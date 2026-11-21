# Fröhlich ballisticity and tunnelling guard, fixed before any device number with them

Written 2026-11-16 (late), before any corrected I_ON was computed. Go-ahead: 16 Nov, after the
Laplace C_G.
- **What it is:** the last two T1 corrections of plan D6: the ballisticity B from a Fröhlich mean free path (LM8,
  LK11), and the tunnelling guard (heavy-plan M16, carried as "WKB guard").
- **Status:** labelled corrections and flags beside the ballistic numbers. Nothing here changes the ballistic T1 or
  Laplace results; the HfS2 benchmark never uses this.
- **Freezing:** `fits/PREREG_corr.md5`; code md5s in `fits/PREREG_corr.code.md5`. Any later change is a dated
  amendment below.

## Common inputs
- **Ribbons:** every T1 width ([100] bare, [110] bare, [110] H) and the four calibration DFT ribbons.
- **Bands:** the cached HSE building-block bands (`fits/t1/cache/`, the same 2000-k grid as T1); the DFT points from
  their HSE device Hamiltonians.
- **Masses:** the HSE band-edge masses already stored in `fits/t1/rows` and `fits/t1/dft.json` (|m_CB| for n, |m_VB|
  for p). One parabolic band per carrier.
- **Device:** L = 15 nm, V_DD 0.7 V, I_OFF 10 nA/μm, 300 K, as signed off. "No claims below L = 15 nm."
- **OFF state:** from the ToB itself (parallel-plate C_G, EOT 0.6 nm). The gate voltage where |I| = I_OFF is found by
  root-finding, not by the 381-point sweep. The barrier above the source Fermi level, E_bF = E_edge + U_OFF − μ_S, and
  the source degeneracy δ_S = μ_S − E_edge (μ_S from the doping) are read from it.
  - **Check:** this pin must reproduce the stored ballistic I_ON within 0.5 %.
  - Under ideal electrostatics, E_bF does not depend on the C_G or EOT, so one value per ribbon and carrier serves
    all three geometries.

## 1. Ballisticity B (code: `scnfet/frohlich.py`)
- **Constants:** ε∞ 8.3, ε0 29.5, ħω_LO 84.9 meV (Grümbel 2024, measured).
- **Main: Sohier's 2D Fröhlich model in vacuum (isotropic, PRB 94, 085415, Eqs. 8, 12–14).**
  - Coupling: |g(q)|² A = K_A / (ε12 + r_eff q)², with ε12 = 1 and r_eff = ε∞ t / 2. The sheet thickness is
    t = 4.50 Å: two planes at the bulk (001) spacing.
  - K_A is fixed by Sohier's large-q limit: there the 2D coupling equals the bulk one, which the measured constants
    give. So no Born charges or phonon eigenvectors are needed.
  - Rate: their Eq. 14 (absorption + emission) for a parabolic band.
  - Vacuum is the least environmental screening, so it is the conservative 2D value. A high-k gate raises ε12 and
    lengthens λ. Remote surface-optical phonons of the HfO2 shorten it; they are left out and stated as a limitation.
- **Other end: the bulk 3D Fröhlich rate** (Lundstrom), same constants.
- **From rate to B:**
  - **Mean free path:** every event is treated as randomising. Half of them reverse the velocity in the 1D channel,
    so λ(E) = 2 v(E) τ(E).
  - **Low-energy (kT-layer) value, as the referee item asks:** λ0 = ∫ λ(E) e^{−E/kT} dE / kT, the one-mode 1D
    Landauer average.
  - **Ballisticity:** B = λ0 / (λ0 + L), the signed-off form. **I_ON,B = B × I_ON** for parallel plate, single gate
    and GAA, at both EOTs, with the IRDS verdict at EOT 0.6.
- **Flags:**
  - **α_3D = (1/ε∞ − 1/ε0) √(m Ry/ħω) > 1:** intermediate coupling, so the perturbative B is not quantitative. This
    applies to the heavy edge bands.
  - **W < 2 nm:** the sheet rate is used without the ribbon's form factor.
- **In-sample checks before freezing** (generic masses, no device numbers):

  | | λ0 at m = 0.2 | m = 0.3 | m = 0.6 | m = 1.0 |
  |---|---|---|---|---|
  | 2D vacuum (main) | 24.5 nm | 16.8 | 9.2 | 6.1 |
  | 3D bulk | 74.6 | 49.8 | 24.9 | 14.9 |
  | 1D backscattering only, form factor 1 (cross-check) | 22.1 | 16.8 | 10.7 | 7.8 |

  - **The bulk 3D rate gives the longer λ.** So it is the optimistic end, not the conservative one the planning
    evidence expected ("bulk formula is expected to give a conservative (short) λ … confirm numerically"). Confirmed
    numerically the other way round.
  - **The strict 1D backscattering estimate agrees with the main one** within 25 %. So the sheet rate with the
    randomising rule is a fair stand-in for the ribbon.
  - **Tests:** the 3D code gives 113 fs with Mu's LDA inputs (Mu ab initio: 104–114 fs). The 2D angular integral
    matches a brute-force q-grid within 3 %. The large-q identity holds exactly.

## 2. Tunnelling guard (code: `scnfet/guard.py`)
- **Model:**
  - effective mass, one band, electron-like (holes mirrored);
  - natural-length barrier ψ'' = (ψ − U_g)/λ_s² on [0, L] with flat leads, max ψ = δ_S + E_bF;
  - λ_s = 1, 2 and 3 nm (M16's range). λ_s = 2 nm is the decision value.
- **Source–drain tunnelling (SDT):**
  - R = I(T_quantum) / I(T_classical) in the OFF state, with Landauer occupations;
  - T_quantum from a transfer matrix (exact in the model). Where T_WKB < e⁻³⁰ the WKB value is used, because the
    matrix loses precision there. The WKB R is reported alongside;
  - **flag if R > 1.5 at λ_s = 2 nm** (M16);
  - ΔV_T = (kT/q) ln R.
  - **Re-pinned I_ON at λ_s = 2 nm:** the thermionic I_OFF target is lowered to I_OFF / R, so that the total OFF
    current is 10 nA/μm. This is done for all three C_G at EOT 0.6, for every ribbon.
  - **Fully corrected I_ON = B × re-pinned I_ON**, with the IRDS verdict.
- **Band-to-band tunnelling:**
  - **Elastic through-path:** in the ideal ToB device, source band → drain opposite band needs E_g < qV_DD − δ_S.
    The margin E_g − qV_DD is reported for every ribbon.
  - **Drain-side window (GIDL-type):** in the OFF state the channel's opposite band, raised by the gate, lines up
    with empty drain states when ΔE = E_bF + qV_DD − E_g > 0.
  - **Upper bound on that current:** two-band Kane WKB through the drain-side profile. It uses m_r = m_c m_v /
    (m_c + m_v) from the band-edge masses, assumes the opposite band is filled, and counts every generated carrier
    as leaving. A floating ribbon cannot remove them, so the true current is lower. The floating-body effect itself
    is outside the model.
  - **Flag if the bound exceeds 0.5 × I_OFF·W at λ_s = 2 nm**, the same 1.5× rule as SDT.
  - **The registered gap criterion (E_g ≥ 0.9 eV, signed off 11 Nov) is reported as it stands.**
- **In-sample checks before freezing** (generic inputs, E_bF 0.30 eV, δ_S 0.05 eV):
  - **SDT R by transfer matrix** (WKB in brackets):

    | | λ_s = 1 nm | λ_s = 2 nm | λ_s = 3 nm |
    |---|---|---|---|
    | m = 0.1 | 0.82 (1.30) | 1.45 (2.07) | 2.20 (2.89) |
    | m = 0.2 | 0.90 (1.17) | 1.21 (1.55) | 1.43 (1.84) |
    | m = 1.0 | 0.99 (1.07) | 1.04 (1.19) | 1.08 (1.26) |

    R < 1 at λ_s = 1 nm is reflection above a sharp barrier. WKB has no such reflection, so it runs high.
  - **Drain-side bound** (m_r 0.2, E_bF 0.29 eV, λ_s 2 nm): 4.3e-9 A at E_g 0.90 eV, against I_OFF·W = 3e-11 A for
    3 nm. It falls to the thermal tail (5.7e-11 A) at E_g 1.0 eV, where the window closes.
  - **Tests:** the square barrier is exact; the heavy-mass R → 1.

## Outputs
- **Data:** `fits/t1/corr/` (one JSON per ribbon, plus `summary.md` and `README.md`), from `scripts/t1_corr.py`.
- **Figure 7 gains the B-corrected curves** (main B; the 3D end as a band). The guard flags go in the table and
  README.
- **Report everything; nothing is tuned.**
