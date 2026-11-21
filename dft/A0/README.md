# A0 — cutoff and pseudopotential tests

Plan §4, stage A0. The criteria below were written down on 2026-11-10, **before any output existed**. The
cutoffs are chosen by these criteria, not by looking at the results.

## Pseudopotentials (SSSP 1.3.0 PBE, md5-verified against the official JSON)
All four files live in `/path/to/scn_fet/pseudo/`. The tarballs and JSON are in
`_sssp_1.3.0_download/` (archive.materialscloud.org record rcyfm-68h65).

| Role | Sc | N | H |
|---|---|---|---|
| PBE relaxations (SSSP precision) | `Sc.pbe-spn-kjpaw_psl.0.2.3.UPF` (PAW, 11 e⁻, 90/720) | `N.oncvpsp.upf` (Dojo NC, 5 e⁻, 80/320) | `H_ONCV_PBE-1.0.oncvpsp.upf` (SG15) |
| HSE06 + W90 (all-NC) | `Sc_ONCV_PBE-1.0.oncvpsp.upf` (SG15 NC, 11 e⁻) | `N.oncvpsp.upf` | `H_ONCV_PBE-1.0.oncvpsp.upf` |

Both Sc files carry 3s3p semicore, so the PAW and NC sets are semicore-consistent.

## Test matrix (`make_a0.py`, 75 runs, one 64-rank job)
| Block | Sets | Runs |
|---|---|---|
| Bulk rocksalt E(V), a = 4.42–4.60 Å (7 points), 12³ Γ-centred, m-v 0.01 Ry | paw60/75/90/110 (dual 8), nc80/90/110 (dual 4) | 49 |
| Bulk at fixed a = 4.50 Å | all 7 | 7 |
| Bulk k ladder at a = 4.50 Å, k = 8/16/20 | paw90, nc80 | 6 |
| Sheet M-ScN, fixed geometry (a 3.11, interlayer 2.18, c 22.5 Å), 12×12×1 | all 7 | 7 |
| Sheet bands Γ–X̄–M̄–Γ | paw90, nc80 | 2 |
| Sheet 2Dxy vc-relax | paw60/75/90/110 | 4 |

nc110 is the reference for the NC set. It is added beyond the plan's 80/90, because a convergence test needs a
point above the candidates.

## Pre-registered criteria
A cutoff **passes** when, measured against the highest cutoff of the same set:
- the bulk lattice constant (Birch–Murnaghan fit) is within 0.002 Å;
- E(V) differences are within 1 meV/atom across the 7 points (curve shape, not absolute E);
- the bulk Γ–X gap at a = 4.50 Å and the sheet gap at fixed geometry are each within 10 meV;
- for PAW only, the sheet a from the 2Dxy vc-relax is within 0.003 Å;
- the pressure at a = 4.50 Å is within 1 kbar.

**Production choice:** the lowest passing cutoff.
- For PAW, a choice below SSSP's 90/720 needs every criterion met.
- For NC, 80 Ry must pass for the plan's HSE setting to stand. If it fails, the HSE cutoff goes up and the
  budget is re-costed (risk register).

**Physics checks.** These are flagged if missed, never tuned to:
- bulk PBE a ≈ 4.50–4.52 Å;
- bulk PBE Γ–X ≈ −0.04 eV (sign and size);
- sheet PBE gap ≈ 0.388 eV (Alexandria);
- PAW vs NC sheet gap at the same geometry within 30 meV;
- the k ladder: gaps converge to within 5 meV by k = 12.

## Run
```bash
python3 make_a0.py                      # local: writes inputs/ and tags.txt
# copy dft/A0/ to /path/to/scn_fet/A0/, then on the cluster:
mkdir -p outputs && sbatch run_a0.sh    # restartable: finished tags are skipped
```
