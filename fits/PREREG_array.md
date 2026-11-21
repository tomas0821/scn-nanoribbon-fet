# Pre-registration: gate capacitance of a ribbon array (closes the Figure 8 bracket)

Written 2026-11-17, before any device number with an array C_G. Go-ahead: 17 Nov, on the optional
item in `fits/t1/pitch/README.md`. Later changes are dated amendments below; this text is not edited.

## What it is
Figure 8 brackets the per-footprint I_ON between parallel plate and the isolated ribbon. This computes the C_G of
the actual array at the two registered spacings (`fits/PREREG_t1b.md` §1), so the bracket collapses to one value per
spacing.

## Geometry (solver: `scnfet/cg_array.py`, using `scnfet/cg.py`'s finite-volume solver and grid rules)
- **The array:** infinitely many identical ribbons at the same potential, with pitch P = W + s.
  - s = 2.0 nm and s = 7.69 nm (= 2 t_HfO2 at EOT 0.6), as registered for Figure 8.
- **The stack:** exactly the `PREREG_cg.md` single gate, so only the neighbours change:
  - a planar gate, with the ribbon's bottom face at t = EOT × 25/3.9 = 3.85 nm;
  - encapsulated in HfO2, ε 25;
  - d = 0.217 nm.
- **One period with zero-flux planes at x = ±P/2** is exact by symmetry. The top boundary is at 40 × max(P, t).
- **GAA needs no new solve.** For s ≥ 2t each ribbon's gate box encloses it, so the array value is the isolated GAA
  value already used. The GAA column of Figure 8 at s = 7.69 nm is therefore already exact, for a zero-thickness
  metal wall between ribbons.

## Numerics (in-sample checks, C_G only, before freezing)
- **s = 0:** C / C_pp = 1.000 for W = 1.09, 2.41 and 5.93 nm. A continuous sheet is the parallel-plate value
  exactly.
- **s = 100 nm:** within 0.3–0.5 % of `cg.cg_single`, the isolated ribbon.
- **Monotonic in s** over 0, 0.5, 1, 2, 4, 7.69, 16, 32 and 100 nm, at all three widths.
- **Convergence:** h 0.02 → 0.01 nm changes C by ≤ 0.024 %; doubling the top extent by ≤ 0.0001 %.
- **What the checks already show:** at s = 2 nm, C / C_pp is 2.66 / 1.75 / 1.31 at those widths (isolated
  7.23 / 4.14 / 2.46). The array sits much nearer parallel plate. This is a C_G number, not a device number.

## Device numbers
- **The same as the Laplace run (`PREREG_cg.md`):**
  - the ToB code, conventions, I_OFF pin and doping;
  - the cached HSE building-block bands (`fits/t1/cache/`);
  - only C_Σ = C_G changes, per unit length per ribbon.
- **Scope:** EOT 0.6 only (Figure 8 and IRDS are at EOT 0.6); n and p; every T1 width of [100] bare, [110] bare and
  [110] H; the four calibration DFT points.
- **B:** × B from `fits/t1/corr/` (2D main, 3D for the range). B does not depend on C_G, so it is reused, not
  recomputed. The source–drain-tunnelling re-pin (≤ 2 %) is not redone.
- **Reported:**
  - C_array / C_pp and C_array / C_isolated;
  - I_ON per W and per footprint (× W/P), ballistic and × B;
  - IRDS pass counts.
- **Figure 8:** the array line goes into the two per-footprint columns. Parallel plate and isolated stay as the
  bracket.
- **The headline verdict stays per W, as registered.** Report everything; nothing is tuned.

## Code and resources
- **`scripts/t1_cg_array.py`** (modes `table`, `tob S` with S = s2 | sGAA) → `fits/t1/cg_array/`.
- **Code md5s** go in `fits/PREREG_array.code.md5` before the first device run.
- **Laptop:** two ToB processes, one per spacing, each with 2 threads and ≤ 2.5 GB virtual memory (the Laplace run
  peaked at 0.29 GB).

## Amendments
*(none)*
