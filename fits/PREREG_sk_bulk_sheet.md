# Pre-registered: the bulk + sheet SK fit and gate (a)

Written 2026-11-11, **before any A2 sheet output existed** (job 141971 was running).

Decisions behind it (2026-11-11):
- **SK form: the plan's 13 parameters plus the three-centre Sc–N–Sc term** (`sk.ThreeCentre`; t3_s, t3_ppp,
  t3_psp; t3_pss fixed at 0). No Sc 4s.
- A2 geometry: scaled by the bulk ratio 4.501/4.5134.

Evidence (2026-11-11): the SK projection and first bulk fit, and the LR1 tests, in `fits/A1_bulk/`.

## Model
- **Orbitals and shells:** Sc d + N p; NN pd, 2NN dd/pp, 3NN Sc–Sc dd, plus the three-centre term on 2NN Sc–Sc
  through every N bonded to both.
- **Harrison exponents:** default `sk.DEFAULT_ETA`; 3.5 per three-centre leg. These are not fitted here: fitting
  them to deformation potentials is a separate, later item.
- **Sheet on-site blocks** (plan D3 "surface on-site blocks"): orbital-resolved shifts by the C4v irreducible
  representations of the sheet sites, the same for every sheet atom:
  - Sc: z² (a1), x²−y² (b1), xy (b2), xz/yz (e);
  - N: pz (a1), px/py (e).

  That makes 6 extra parameters. They apply only to the sheet structure; bulk is unchanged.
- **Energy gauge:** Ep = 0 in bulk. The sheet is referenced to its own VBM, through its own `ref`.

## Targets and weights (fixed now)
**Bulk:** exactly as `scripts/a1_sk_fit.py`:
- HSE06 W90 16³ path ≤ VBM + 4 eV, plus the DFT 16³ grid eigenvalues up to VBM + 6.5 eV, with floors above;
- characters X3 = d_xy and Γ25′ = t2g.

**Sheet:** the same recipe on the A2 data:
- W90 16×16 path Γ–X̄–M̄–Γ ≤ VBM + 4 eV;
- DFT 16×16 IBZ eigenvalues of the 16 model bands up to the frozen-window top (A2-2), with floors above.

**Weights (the "Tn" recipe):**
- σ 20 meV near the gap (VB ≥ −1.5 eV; CB ≤ CBM + 2.5 eV), 0.5 eV elsewhere;
- masses at relative σ 0.01 (electron) and 0.0167 (holes), i.e. ×3.

**Masses:**
- bulk: m_l, m_t at X; heavy and light hole along [100] and [110];
- sheet: the conduction and valence band-edge masses at the band-edge k along Γ–X̄ and Γ–M̄ (both bands at each
  edge if degenerate).

**Seeds:** the bulk Tn fit, with the sheet on-site shifts at 0, plus the projection of the sheet Wannier H(R).
Multi-start as before.

## Gate (a), as plan §4, made explicit
Judged on this one joint fit, with no refit after seeing the score:
- RMS within ±1 eV of the band edges (VB ≥ VBM − 1 eV; CB ≤ CBM + 1 eV) ≤ **30 meV**, separately for bulk and sheet
  (the worse of path and grid each);
- **every transport mass listed above within ±10 %** of its Wannier value.

**Pass:** proceed to A3 and the edge blocks.

**Fail:** report the failure. Plan §10 LR1 then applies:
- the remaining pre-registered mitigation, Sc 4s, was already tested on bulk and did not help the masses;
- so the next step is a decision between DFT-Wannier device Hamiltonians for 3–4 ribbons (+3–5k core-h) and
  carrying the SK mass errors into the device error bars;
- weights will not be tuned to make the gate pass.

The PBE twin (D7) uses the same form, recipe and targets from the PBE runs.
