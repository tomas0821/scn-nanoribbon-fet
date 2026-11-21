# Building-block Wannier ribbons on the PBE hold-outs (plan LR2; laptop, 2026-11-16)

Recipe and criteria: `fits/PREREG_bb_holdouts.md`, frozen at 16:20 on 16 Nov (md5 `8e722437…`, checked by the run;
code md5s in `fits/PREREG_bb_holdouts.code.md5`), before any building-block Hamiltonian of a hold-out was assembled.

- **Code:** `scnfet/bb.py` (the assembler, with `scnfet/tests/test_bb.py`; 58 tests pass) and
  `scripts/gate_bb_holdouts.py` (modes `selftest`, `calib`, `run`).
- **Outputs:** `selftest.json`, `calib.json`, `gate_bb.json`, `run.log`, and `figures/Fig3b_holdouts_building_blocks.png`.
- **Resources:** capped at 2 threads and 2.5 GB, with peak RSS 475 MB. The run took 12 min.

**Recipe in one line:** elements between atoms within the calibration ribbon's half-width of the same edge come from
the calibration ribbon's Wannier H (variant W1 / HW1, the sheet's window rule); every other element is the sheet's,
rotated into the ribbon's axes; the sheet's on-site energies are matched to the calibration ribbon's centre rows;
the cutoff is 20 Å; nothing is fitted.

## Result: **both hold-outs pass every registered item**

**`r110_8` ([110] n = 8, W 2.31 nm; edges from `r110_4`)**

| | PBE-Wannier W2 | building blocks | Δ | limit | |
|---|---|---|---|---|---|
| gap | 0.285 eV | 0.252 eV | − 32 meV | ± 50 | ok |
| valley E_c(X) − E_c(Γ) | 1.155 | 1.148 | − 7 meV | ± 50 | ok |
| m_VB / m_CB | − 0.289 / 0.959 | − 0.266 / 0.903 | − 8.0 % / − 5.9 % | ± 10 % | ok |
| I_ON n EOT 1 / 0.6 (mA/μm) | 0.832 / 1.402 | | + 0.9 % / + 0.8 % | ± 15 % | ok |
| I_ON p EOT 1 / 0.6 | 1.449 / 2.317 | | + 3.7 % / + 3.9 % | ± 15 % | ok |
| IRDS verdict | n no, p yes | n no, p yes | same | | ok |
| **HSE gap** (HW1 + HSE sheet) | HSE 0.921 eV | 0.911 eV | **− 10 meV** | ± 80 | **ok** |

**`r100_6` ([100] n = 6, W 2.41 nm; edges from `r100_3`)**

| | PBE-Wannier W2 | building blocks | Δ | limit | |
|---|---|---|---|---|---|
| gap | 0.617 eV | 0.609 eV | − 8 meV | ± 50 | ok |
| valley | 1.637 | 1.635 | − 3 meV | ± 50 | ok |
| m_VB / m_CB | − 0.199 / 0.203 | − 0.198 / 0.204 | − 0.4 % / + 0.5 % | ± 10 % | ok |
| I_ON n EOT 1 / 0.6 | 2.149 / 3.752 | | − 0.3 % / − 0.3 % | ± 15 % | ok |
| I_ON p EOT 1 / 0.6 | 2.082 / 3.426 | | − 1.0 % / − 1.3 % | ± 15 % | ok |
| IRDS verdict | n yes, p yes | n yes, p yes | same | | ok |
| **HSE gap** (HW1 + HSE sheet) | HSE 1.629 eV | 1.605 eV | **− 25 meV** | ± 80 | **ok** |

**Reading the [110] result:**
- The [110] margins are moderate: gap − 32 of ± 50 meV, hole mass − 8.0 of ± 10 %. [100] is within 1 %.
- **The in-sample prediction was wrong in the good direction.** Rebuilding `r110_4` from its own halves plus the sheet
  had moved its hole mass by − 40 % (PBE) and + 88 % (HSE). The n = 8 interior is evidently far more sheet-like than
  the n = 4 centre, which sits only 3–4 rows from polar edges.
- The Sc-s edge conduction band transfers with the edge region: it is the CB here, as in the reference.
- **Without the edge-s bands** the building blocks give a gap of 0.748 eV and a valley of 0.652 eV, against W2's
  0.807 and 0.633. That is descriptive, not registered.
- **The building blocks carry the vacuum-like state.** In [110] they show extra bands at Γ about 1.1 eV above the CBM
  (figure, left). These are the sheet's vacuum-reaching state: W1's window freezes it, while the W2 reference
  deliberately leaves it out. A5-2 put the hold-out's state at CBM + 1.105 eV, so it lies outside the device window.

## Descriptive (no pass/fail)
| | `r110_8` | `r100_6` |
|---|---|---|
| vacuum alignment only (sheet shift 0) | gap + 46 meV, **valley − 107 meV**, p I_ON − 20 %, p IRDS verdict differs → would fail | gap + 2, valley + 20, I_ON ≤ 6 % → would pass |
| W2 edges (gauge-mismatched), centre-aligned | **gap − 109 meV**, m_VB − 16 %, p I_ON + 19–20 % → would fail | gap + 4, valley + 26, masses ≤ 4 % → would pass |

- **[100] passes with any of the three variants. [110] passes only with the registered recipe.**
- Both [110] choices were made before the run, on in-sample grounds:
  - the shared Wannier gauge (W1 rather than W2): with W2, the Sc z² elements disagree with the sheet's by 0.2–0.9 eV;
  - the centre-row alignment: the vacuum-referred narrow-ribbon centre sits 0.1–0.2 eV off the sheet.
- The gate (d) drift result (≈ + 40 meV, [110]) was known when the alignment was chosen (PREREG, "Known before
  freezing"). Report it that way.
- The registered sheet shifts are − 0.212 eV ([100] PBE), − 0.243 eV ([100] HSE), + 0.139 eV ([110] PBE) and
  + 0.196 eV ([110] HSE).

## What it means (by the pre-registered consequence, the 16 Nov decision read per direction)
- **Both directions pass → building-block Wannier ribbons are the T1 device Hamiltonians for bare [100] and bare
  [110] at 1–6 nm.** The fallback is not triggered.
- **H-terminated ribbons:** building-block curves from the H calibration ribbons, labelled "validated on bare edges
  only". The H recipe's own in-sample checks go first.
- **SK stays in the paper** as the fitted model with measured transfer errors: gate (d), `fits/holdouts/README.md`.
  It is overlaid on [100] as a cross-check.
- **For the HSE device Hamiltonians:** the building blocks reproduce both HSE gap points within 25 meV, which is
  the only HSE check available on the hold-outs.
