# T1: building-block device Hamiltonians across 1–6 nm (laptop, 2026-11-16)

Settings: `fits/PREREG_t1.md`, frozen at 16:43 on 16 Nov before any width result (md5 in `.md5`; code md5s in
`fits/PREREG_t1.code.md5`).
- **Code:** `scripts/t1_bb_widths.py` (modes `hcalib`, `geomcheck`, `sweep`, `dft`, `sk`) and `scripts/t1_figs.py`.
- **Outputs:** `rows/<family>_n<n>.json`, `sk/`, `dft.json`, `hcalib.json`, `geomcheck.json`, `t1_table.md`; figures
  `figures/Fig4_Eg_W.png` (with error bands since 17 Nov), `Fig5_valley_mass_W.png` and
  `Fig7_Ion_W_parallel_plate.png`. Figures 7, 7b and 8 come from `cg/`, `corr/` and `pitch/`.
- **Resources:** two background processes, each capped at 2 threads and 2.5 GB; peak RSS 2.1 GB each, in the largest
  eigen-solves. The SK overlay was paused while the bare sweep ran, so the two peaks never overlapped. Wall time
  ≈ 2.5 h.

**All device numbers here are "ballistic, ideal-electrostatics upper bound; parallel-plate C_G".** SS is 59.5 mV/dec
everywhere, as ideal electrostatics gives.
- **Laplace C_G: done 16 Nov, in `cg/README.md`.** With it, every family clears IRDS for n and p (Figure 7).
- **Decided (16 Nov):** parallel plate carries the design rule; the Laplace result is the other bound.
- **Fröhlich B and the tunnelling guard: done 16 Nov, in `corr/README.md` (Figure 7b).** The heavy edge bands get
  B ≈ 0.3, which brings the design rule back under the Laplace C_G. The drain-side band-to-band bound flags [110]
  bare from 1.38 nm and [110] H from 2.63 nm.
- **Pitch-normalised I_ON: done 17 Nov, in `pitch/README.md` (Figure 8; `fits/PREREG_t1b.md` §1).** At s = 2 nm
  and 7.7 nm, the family ranking is the same as per W. Per footprint × B, no ribbon passes IRDS under parallel
  plate. With the isolated single gate at s = 2 nm, the design rule holds at every width.
- **Array C_G: done 17 Nov, in `cg_array/README.md` (Figure 8; `fits/PREREG_array.md`).** With the array's own C_G,
  the design rule holds per W at both spacings.
- **Figure 6 (transfer characteristics): `scripts/t1_fig6.py`.**
- **Error bands: done 17 Nov (Figure 4; `fits/PREREG_t1b.md` §2).** Constant per direction: ± 25 meV [100],
  ± 32 meV [110]; [110] H borrows the [110] band.

## Two code fixes before the sweep (sanity checks only; no number changed)
`scnfet/bb.py` was re-hashed after each, and the reason is in `PREREG_t1.code.md5`.
- **Atom-match tolerance:** relaxed H in the H calibration ribbons sits up to 0.7 Å (y) and 1.6 Å (z) from the ideal
  cut. Matching is by ideal row, so only the sanity tolerance changed: < Lx/4 in x, < 2 Å in y.
- **Real-gauge threshold:** `r110_4h` HW1 carries 0.4 meV of imaginary noise. The threshold went from 0.1 to 1 meV;
  the bare ribbons are all below 0.1 meV.

## Checks before the sweep
- **D4 copy-step geometry:** the copied geometry against the relaxed hold-outs gives W 23.088 / 23.115 Å (`r110_8`) and
  24.099 / 24.120 Å (`r100_6`). The largest atom deviation is 0.031 Å.
- **The H rule** (`hcalib.json`): each H calibration ribbon rebuilt from its halves plus the sheet, against itself.
  - **`r100_3h`:** PBE passes (gap − 8 meV, masses ≤ 3 %). **HSE fails on the hole mass, − 12.1 %** (gap − 13 meV,
    valley − 26 meV). **So [100] H gets no curve: DFT point only.**
  - **`r110_4h`:** PBE gap − 36 meV, masses ≤ 8.2 %; HSE gap − 43 meV, masses ≤ 8.4 %, valley ≤ 4 meV. Passes →
    **a labelled curve, "validated on bare edges only; in-sample for H".**
  - **Caveat for any H curve:** the H calibration ribbons' centre rows sit 2.1–3.0 eV below the vacuum-referred sheet
    (bare: ± 0.1–0.2 eV). The H edges carry strong electrostatics.
    - The centre-row alignment absorbs this at the calibration width.
    - The in-sample check cannot test whether the interior-to-edge potential changes with width.

## Results (HSE06 building blocks; full table in `t1_table.md`)
| | [100] bare | [110] bare | [110] H (in-sample only) |
|---|---|---|---|
| widths | n 3–14, 1.09–5.93 nm | n 4–20, 1.07–6.04 nm | n 4–20, 1.39–6.36 nm |
| HSE gap | 1.952 → 1.452 eV, direct at Γ | 1.029 → 0.877 eV (CB = Sc-s edge pair) | 1.444 → 0.902 eV, indirect (VBM at X: Sc–H edge band) |
| E_c(X) − E_c(Γ) | + 1.75 → + 2.00 eV | + 1.62 → + 1.54 (edge-s CB); interior bands + 0.23 → + 0.94 | + 0.75 → + 0.96 |
| m_CB / m_VB | 0.27 → 0.19 / − 0.19 → − 0.18 | 0.89 → 0.98 (edge s) / flat at n = 4, − 0.39 → − 0.24 | 0.29 → 0.19 / − 0.85 (flat edge band) |
| I_ON n, EOT 0.6 (mA/μm) | 2.76 – 3.80 (peak at 3.3 nm) | 1.34 – 1.61 | 2.78 – 3.82 |
| I_ON p, EOT 0.6 | 3.07 – 3.67 (peak at 2.0 nm) | 1.53 – 2.41 | 1.43 – 1.82 |
| IRDS (≥ 1.979 at EOT 0.6) | n and p at every width | n never; p from 1.69 nm (15 of 17) | n at every width; p never |

- **DFT points agree at the calibration widths.** At EOT 0.6 nm, I_ON (n / p, mA/μm) for the HSE device
  Hamiltonians against the building blocks:
  - [100] n = 3, HW3 (vacuum-like state included): 2.75 / 3.34 against 2.76 / 3.35;
  - [110] n = 4, HW1: 1.31 / 1.48 against 1.35 / 1.53;
  - [110] H n = 4, HW2: 2.66 / 1.42 against 2.78 / 1.43;
  - [100] H n = 3 (DFT point only): 2.87 / 0.81, with the p-FET held back by its Sc–H edge VB band.
- **SK-HSE cross-check on [100]:**
  - **Gap:** 25–95 meV wider than the building blocks (+ 95 meV at 1.1 nm, + 25 meV at 5.5 nm).
  - **n-FET I_ON:** within ≈ 3 %.
  - **p-FET I_ON:** high by up to + 23 % at 1–2.4 nm, within ≈ 3–8 % from 3.3 nm. This is the SK light-hole
    error (gate (a); − 32 % in-sample HSE at n = 3).
- **PBE twin, bands only:** [110] H closes its PBE gap from 2.3 nm upward (− 0.02 eV at n = 7, − 0.20 eV at 6.4 nm). Its
  Sc–H edge VB band overlaps the interior CB as the interior gap shrinks towards the sheet's 0.40 eV. HSE keeps
  0.90 eV, and HSE is the production functional.

## Flags for the claims and the next steps
1. **No conduction-valley reordering anywhere in 1–6 nm (HSE or PBE).** Every E_c(X) − E_c(Γ) stays positive.
   - In [110], the sheet's M̄-derived valley (the interior bands, edge s left out) closes from + 0.94 eV at 6 nm to
     + 0.23 eV at 1.07 nm. The DFT calibration ribbon itself (n = 4) shows the same, so this is a DFT result within
     the bracket.
   - A crossover would need W < 1.1 nm. That is outside both the range and the building blocks (n ≥ 4), so the
     paper cannot claim one, even as a model prediction from these Hamiltonians.
   - **This changes a headline result of the plan.** **Decided (16 Nov):** write "the [110] M̄ valley
     approaches to within 0.23 eV of Γ at 1.1 nm without crossing" as a supporting result. The headline becomes
     the orientation/termination design rule of flag 2.
   - In bare [110], the lowest CB is the Sc-s edge pair anyway (LR3).
2. **[110] bare n-FET fails IRDS at every width.** The heavy Sc-s edge CB (m ≈ 0.9–1.0) gives 1.3–1.6 mA/μm. H
   termination removes the edge CB and puts the n-FET at [100] level, but its Sc–H edge VB sinks the p-FET.
   - **Holds under parallel plate only.** With the Laplace C_G, both edge-band penalties shrink (to 0.6–0.87 of
     [100]), and every family passes (`cg/README.md`).
3. **BTBT guard:** the plan's 0.9 eV criterion. The [110] bare HSE gap drops below 0.9 eV from 2.9 nm (0.877 eV at
   6 nm). [110] H stays at ≥ 0.902 eV. The WKB guard decides.
   - **Guard result (`corr/README.md`):** there is no elastic band-to-band path (all gaps exceed qV_DD = 0.7 eV). The
     drain-side window bound flags [110] bare from 1.38 nm and [110] H from 2.63 nm.
4. **The [100] I_ON(W) hump:** n peaks at 3.3 nm and p at 2.0 nm, both inside ballistic ToB. With the Laplace C_G
   it is gone: I_ON falls monotonically with W.
5. **T1 is complete (17 Nov).** *(Laplace C_G done 16 Nov: `cg/README.md`, Figure 7. Fröhlich B and the guard
   done 16 Nov: `corr/README.md`, Figure 7b. Pitch-normalised I_ON and the error bands done 17 Nov: `pitch/README.md`,
   Figures 8 and 4.)*
   - **The thickness axis is dropped** (decision, 16 Nov). The 3-plane slab stays as a DFT data point.
   - Gate (e), the HfS2 benchmark, needs QUDEN.

## Post-hoc diagnostic (18 Nov, not registered; no verdict changes): where the H self-rebuild deviates
`scripts/diag_hcalib_dev.py` → `hcalib_dev.json`. Same code path as `hcalib`. Both H ribbons have an even row count,
so the rebuild copies every on-site and same-half element from the ribbon itself. Only the couplings across the
centre come from the sheet, and the sheet shift never enters.

| HSE06 | ±1 eV max (sorted index) | device window (±0.5 eV) max | VB within 0.6 eV of VBM | I_ON rebuild / own (EOT 0.6) |
|---|---|---|---|---|
| [110] H | 312 meV, CB+1 at CBM + 0.99 eV, k = 0.081 (99 % centre rows, 92 % Sc d) | 43 meV, at the CBM (= the gap error) | ≤ 0.5 meV (the Sc–H edge band) | n +4.9 %, p +0.2 % |
| [100] H | 285 meV, CB+2 at CBM + 0.89 eV, k = 0.177 (95 % centre rows, 89 % Sc d) | 28 meV (VB−3 at Γ, 0.2–0.4 eV below the VBM) | ≤ 28 meV | n −0.2 %, p +4.5 % |

- **[110] H conduction side:** the deviation falls from 43 meV at the CBM to 18 meV at CBM + 0.4–0.5 eV. It jumps to
  150–312 meV from CBM + 0.5 eV, in the higher interior Sc d subbands.
- **[100] H conduction side:** ≤ 15 meV up to CBM + 0.6 eV, then 203–285 meV from CBM + 0.7 eV. An overlap-matched
  pairing gives 390 meV on CB+0 at CBM + 0.78 eV (overlap 0.83), also outside the device window.
- **PBE:** the same pattern, smaller. [110] H gives 53 meV on CB+1 at CBM + 0.77 eV; [100] H gives 62 meV on CB+3 at
  CBM + 0.98 eV.
- **Reading:** the large numbers are the interior conduction subbands 0.8–1.0 eV above the CBM. They are set by the
  sheet's cross-centre Sc d couplings. They are not edge bands, and they lie outside the device window.
