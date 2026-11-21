# Claims gates (b), (c), (d): the SK models on the hold-outs (laptop, 2026-11-16)

Recipe and criteria: `fits/PREREG_sk_holdouts.md`, frozen at 13:04 on 16 Nov (md5 `c2f5eb0a…`, checked by the run)
before any SK hold-out number existed.

Script: `scripts/gate_bcd_sk_holdouts.py` with modes `selftest`, `calib`, `run` and `figs`. Outputs: `selftest.json`,
`gate_c.json`, `gate_bd.json`, `run.log`, and figures `figures/Fig3_holdouts_PBE_twin.png` and
`figures/Fig2c_slab_holdout.png`.

Runs were capped at 2 threads and 2.5 GB, with peak RSS ≈ 230 MB. The full run takes about 15 minutes, most of it
in the ToB.

**One numerical fix between the frozen recipe and the result.** The first `run` stopped at the first Wannier I_ON
call (`run_attempt1_sweep.log`), before any SK number was computed:
- ToB's gate sweep of ±1.5 V did not bracket I_OFF. The registered doping (0.01 per atom, 4.1e13 cm⁻² in M-ScN) needs
  q·n_SD/C_G ≈ 2 V to deplete at EOT 1 nm.
- `ToBFET.metrics` now takes `span`. The default is unchanged, and a new test shows a wider sweep gives the same I_ON
  and I_OFF (55 tests pass).
- The run uses ±4 V at the same ≈ 21 mV step. The sweep only brackets the I_OFF pin, so no convention changed.

## Gate (b): 3-plane slab — **PBE twin FAILS; SK-HSE gap passes**
| | PBE-Wannier W2 | PBE-twin SK | Δ | limit | |
|---|---|---|---|---|---|
| gap | 0.120 eV | 0.175 eV | + 54 meV | ± 50 | fail (just) |
| E_c(M̄) − E_c(Γ) | 0.622 | 0.620 | − 2 meV | ± 50 | ok |
| E_c(X̄) − E_c(Γ) | 3.630 (DFT 3.143) | 1.911 | − 1,720 meV (vs DFT − 1,232) | ± 50 | fail |
| m_VB Γ–X̄ / Γ–M̄ | − 0.354 / − 0.188 | − 0.390 / − 0.397 | + 10 % / + 112 % | ± 10 % | fail |
| m_CB Γ–X̄ / Γ–M̄ | 0.089 / 0.086 | 0.155 / 0.152 | + 73 % / + 77 % | ± 10 % | fail |
| **SK-HSE gap** | HSE grid 1.1055 eV | 1.136 eV | **+ 30 meV** | ± 80 | **ok** |

- The two top valence bands are degenerate at Γ, so the band-index VB masses mix them. VB−1: + 3 % / − 25 %.
- The SK X̄ failure is a spurious flat conduction band near 2 eV along X̄–M̄ (figure, left panel). The middle plane
  has bulk on-site energies but sheet-like hoppings, so this is the thickness-axis transfer failing in the e_g
  manifold, as the gate (a) analysis predicted.
- The M̄ offset, the quantity behind the valley-crossover claim (LR4), transfers to 2 meV.
- Masses here use one fit rule for both models, on the 601-point path. W2's CB mass is 0.089, against 0.099 in
  `a3_criteria.json`, which used the DFT path's fit window.

## Gate (d): hold-out ribbons
**`r100_6` ([100] n = 6, W 2.41 nm): FAILS on the hole mass only.**

| | PBE-Wannier W2 | PBE-twin SK + edge blocks | Δ | limit | |
|---|---|---|---|---|---|
| gap | 0.617 eV | 0.624 eV | + 6 meV | ± 50 | ok |
| valley E_c(X) − E_c(Γ) | 1.637 | 1.625 | − 12 meV | ± 50 | ok |
| m_VB | − 0.199 | − 0.162 | **− 18.6 %** | ± 10 % | **fail** |
| m_CB | 0.203 | 0.223 | + 9.7 % | ± 10 % | ok (at the limit) |
| I_ON n, EOT 1 / 0.6 (mA/μm) | 2.149 / 3.752 | 2.045 / 3.596 | − 5 % / − 4 % | ± 15 % | ok |
| I_ON p, EOT 1 / 0.6 | 2.082 / 3.426 | 2.207 / 3.474 | + 6 % / + 1 % | ± 15 % | ok |
| IRDS verdict (≥ 1.979 at EOT 0.6) | n yes, p yes | n yes, p yes | same | | ok |
| **SK-HSE gap** | HSE 1.629 eV | 1.686 eV | **+ 57 meV** | ± 80 | **ok** |

**`r110_8` ([110] n = 8, W 2.31 nm): FAILS on every item. This is LR3, as gate (c) predicted.**

| | PBE-Wannier W2 | PBE-twin SK + edge blocks | Δ | |
|---|---|---|---|---|
| gap | 0.285 eV | 0.782 eV | + 497 meV | fail |
| valley E_c(X) − E_c(Γ) | 1.155 | 0.700 | − 455 meV | fail |
| m_VB / m_CB | − 0.289 / 0.959 | − 0.177 / 0.279 | − 39 % / − 71 % | fail |
| I_ON n EOT 1 / 0.6; p EOT 1 / 0.6 | 0.832 / 1.402; 1.449 / 2.317 | 1.776 / 3.157; 2.094 / 3.333 | + 113 / + 125 %; + 45 / + 44 % | fail |
| IRDS n | no | yes | differs | fail |
| **SK-HSE gap** | HSE 0.921 eV | 1.777 eV | **+ 856 meV** | **fail** |

- The Wannier CB is the Sc-s edge pair, which SK has no orbital for.
- **Without the edge-s bands**, W2's gap is 0.807 eV and its valley offset 0.633 eV. SK's interior conduction band is
  then 25 meV low in gap terms and + 67 meV in valley terms (descriptive; not a registered comparison).

**Without the edge blocks (descriptive):**
- `r100_6`: gap + 21 meV, valley − 9 meV, masses − 20 % / + 11 %, I_ON within 8 %.
- `r110_8`: gap + 528 meV.
- So the [100] blocks bring the gap from + 21 to + 6 meV and change little else.

**Edge-block transfer (the D4 drift check):** the hold-outs' own W2 edge blocks differ from the calibration ribbons'
by at most 0.048 eV (Sc) and 0.006–0.012 eV (N). The edge on-site blocks do not drift with width (n = 3 → 6 and
4 → 8).

**ToB k-sampling:** 4000 against 2000 k changes the Wannier I_ON by ≤ 1e-3.

## Gate (c): in-sample, before freezing (`gate_c.json`)
SK + its own edge blocks against the calibration Wannier Hamiltonian:

| Ribbon | Functional | gap error | masses |
|---|---|---|---|
| [100] n = 3 | PBE | + 42 meV | CB ≤ 1.1 %; VB − 16 % |
| [100] n = 3 | HSE | + 82 meV | CB ≤ 1.1 %; VB − 32 % |
| [110] n = 4 | PBE | + 1,040 meV | (Sc-s edge CB) |
| [110] n = 4 | HSE | + 1,314 meV | (Sc-s edge CB) |

The in-gap states were characterised in A4a and A4b: the [110] bare Sc-s edge pair (LR3), and the [100] H and
[110] H Sc–H valence edge bands.

## What the gate says (judged by the frozen criteria)
- **(a)** failed (11 Nov), LR1 → DFT-Wannier device Hamiltonians (decided).
- **(b) fails** for the PBE twin: gap by 4 meV over the limit, the X̄ valley, and the masses. The SK-HSE gap passes.
- **(c)** The edge blocks are calibrated and do not drift with width. The bare [110] Sc-s edge state is outside the
  SK basis.
- **(d) fails for both hold-outs.**
  - `r100_6` misses only on the hole mass (− 18.6 %, inherited from gate (a)'s sheet hole masses, − 30 to − 33 % in the PBE twin). Its
    gap, valley, electron mass, I_ON, IRDS verdict and HSE gap all pass.
  - `r110_8` fails throughout, because SK cannot hold the Sc-s edge conduction band.
- **(e)** is still open: the HfS2 ToB benchmark needs the QUDEN download.
- **Per the pre-registration: plan LR2.** No refitting; the failures are reported as transferability limits. Which
  device Hamiltonians T1 uses at each width was decided separately (`fits/bb/README.md`).
