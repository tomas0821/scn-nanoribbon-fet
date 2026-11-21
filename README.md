# scn-nanoribbon-fet

Code and data for the paper

> **Edge orientation and termination in scandium nitride nanoribbon transistors: ballistic limits from Wannier
> building blocks**
> Tomas Rojas and Ruhi Thorat (manuscript in preparation; the DOI will be added on publication).

The repository holds the `scnfet` Python package and the scripts that turn deposited DFT and Wannier90 results into
every laptop-stage number and figure of the paper:

- **building-block Wannier Hamiltonians** for M-ScN nanoribbons along [100] and [110], 1–6 nm wide, bare and
  H-terminated: calibration-ribbon edges plus the periodic sheet's interior, after Lee, Nardelli and Marzari, PRL 95,
  076804 (2005);
- **ballistic top-of-barrier (ToB) n- and p-FET metrics** (I_ON at fixed I_OFF, SS) against IRDS targets;
- three corrections: the **geometric gate capacitance** (Laplace solve: isolated single gate, gate-all-around and a
  periodic ribbon array), a **2D Fröhlich ballisticity** bound and a **tunnelling guard**
  (source–drain WKB and a band-to-band bound);
- a **Slater–Koster (SK) model** of 2D ScN (Sc 3d + N 2p, up to third-neighbour Sc–Sc, with a three-centre Sc–N–Sc
  term) fitted to HSE06 bands of bulk ScN and the M-ScN sheet, and tested on held-out DFT structures;
- the hold-out tests of both models, and the HfS2 benchmark of the ToB code.

M-ScN, the two-plane rocksalt-(001) bilayer, is a **hypothetical** material: two-dimensional ScN has not been
synthesised. All device metrics are ballistic, ideal-electrostatics upper bounds.

## Layout

| Path | Contents |
|---|---|
| `scnfet/` | the package: `w90` (Wannier90 `_hr.dat` I/O, band interpolation), `bb` (building blocks), `tob` (top-of-barrier FET), `cg`, `cg_array` (Laplace gate capacitance), `frohlich`, `guard` (tunnelling), `sk`, `skfit`, `structures`, `targets_lit`; `tests/` (pytest); `data/` (digitised literature band paths used by the T0 scaffold fit) |
| `scripts/` | one script per analysis stage (below) |
| `dft/` | DFT and Wannier90 inputs, input generators, parsers, job scripts, parsed summaries, the Wannier Hamiltonians the scripts read, and the few raw QE outputs they parse (see `dft/README.md`) |
| `fits/` | results: SK fits (`T0_literature`, `A1_bulk`, `A2_joint`), hold-out gates (`holdouts`, `bb`), HfS2 benchmark (`gate_e`), device results (`t1/`) |
| `fits/t1/cache/` | HSE06 building-block eigenvalues of every T1 ribbon on the 2000-k ToB grid (`.npz`); the C_G, correction, array and Fig. 6 scripts start from them. Not in git: a release asset (see "T1 band cache" below) |
| `figures/` | output directory for the figure scripts (empty in the repository) |

## Installation

Python 3.12. The analysis environment used for the paper:

```bash
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt        # numpy 1.26.4, scipy 1.17.1, matplotlib 3.11.2, pytest 9.2.0
```

The two structure renderings (`scripts/render_ribbons.py`, `scripts/render_structures.py`, Fig. 1) need OVITO and
ASE in a separate environment (`requirements-render.txt`: ovito 3.17.0, ase 3.29.0, pillow 13.0.0). Nothing else in
the repository needs them.

## Running

Every script runs from the repository root:

```bash
PYTHONPATH=. python -P scripts/<script>.py <mode> [args]
python -m pytest scnfet/tests                 # unit tests (71, about 10 s)
```

The docstring at the top of each script lists its modes, inputs and outputs. Scripts that recompute a stored result
overwrite it in `fits/`; figures go to `figures/`. Most modes take seconds to minutes on a laptop. The width sweep
(`t1_bb_widths.py sweep`) and the cache build (`t1_cg.py cache`) peak at about 2 GB per process for the widest
ribbons (n ≈ 14–20); run one such process at a time. The paper's runs used 2 BLAS threads per process
(`OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2`).

### 1. SK model (bulk + sheet)
| Step | Command | Output |
|---|---|---|
| scaffold fit on literature bands | `t0_sk_literature_fit.py` | `fits/T0_literature/` |
| Wannier → SK linear projection (bulk) | `a1_sk_projection.py` | `fits/A1_bulk/sk_projection.json` |
| bulk fits and model-extension tests | `a1_sk_fit.py`, `a1_sk_fit_extensions.py`, `a1_sk_fit_s.py`, `a1_sk_fit_t3.py`, `a1_sk_tradeoff_plot.py` | `fits/A1_bulk/` |
| joint bulk + sheet fit (gate a) | `a2_sk_fit_bulk_sheet.py HSE06`, `… PBE` | `fits/A2_joint/sk_joint_{HSE06,PBE}.json` (the SK parameters) |
| diagnostics of the joint fit | `a2_sk_failure_diagnostics.py` | `fits/A2_joint/` |

### 2. Hold-out tests and benchmark
| Step | Command | Output |
|---|---|---|
| SK on the slab and ribbon hold-outs (gates b–d) | `gate_bcd_sk_holdouts.py selftest`, `calib`, `run`, `figs` | `fits/holdouts/` |
| building blocks on the ribbon hold-outs | `gate_bb_holdouts.py selftest`, `calib`, `run`; figure: `fig_holdouts_bb.py` | `fits/bb/`, `figures/Fig3b_holdouts_building_blocks.png` |
| HfS2 ToB benchmark (gate e) | `gate_e_hfs2_benchmark.py` (needs QUDEN_MATERIALS, below); `--plot-only` redraws the figure from the stored result | `fits/gate_e/` |

### 3. Device results (T1)
| Step | Command | Output |
|---|---|---|
| H-termination rule, geometry check | `t1_bb_widths.py hcalib`, `geomcheck` | `fits/t1/hcalib.json`, `geomcheck.json` |
| width sweep (HSE bands + ToB, PBE bands) | `t1_bb_widths.py sweep 100 110 110h` | `fits/t1/rows/` |
| DFT overlay points; SK overlay on [100] | `t1_bb_widths.py dft`, `sk` | `fits/t1/dft.json`, `fits/t1/sk/` |
| Laplace C_G table; band cache; ToB with C_G | `t1_cg.py table`, `cache`, `tob single`, `tob gaa` | `fits/t1/cg/`, `fits/t1/cache/` |
| Fröhlich B and tunnelling guard | `t1_corr.py run` (families and DFT points) | `fits/t1/corr/` |
| pitch-normalised I_ON | `t1_pitch.py` | `fits/t1/pitch/` |
| periodic-array C_G | `t1_cg_array.py table`, `tob s2`, `tob sGAA`, then `t1_array_post.py` | `fits/t1/cg_array/` |
| transfer characteristics | `t1_fig6.py` (`--plot-only`: figure from the stored `fig6.json`) | `figures/Fig6_transfer.png`, `fits/t1/fig6.json` |
| figures (PNG and vector PDF, lettering ≥ 8 pt) | `t1_figs.py figs`: every T1 figure of the paper and SM (`Fig4_Eg_W`, `Fig5_valley_mass_W`, `Fig7_Ion_W`, `Fig7b_Ion_W_corrected`, `Fig8_Ion_pitch`); writes nothing in `fits/` | `figures/` |
| tables | `t1_figs.py`, `t1_figs.py corr`, `t1_figs.py pitch` (older modes: rewrite the `fits/t1/*.md` tables and draw earlier figure versions) | `fits/t1/*.md` |
| diagnostic of the H-rule deviation | `diag_hcalib_dev.py` | `fits/t1/hcalib_dev.json` |
| post hoc checks, not pre-registered (B under oxide screening and the charge-pinned form; Table 3 margins; HSE06 building blocks at the hold-out SCF k-points) | `posthoc_b_variants.py`, `posthoc_ion_margins.py`, `posthoc_hse_holdout_kcheck.py` (needs the A5 HSE06 outputs), `posthoc_pbe_chord.py` | `fits/posthoc_2026-11-19/` |
| post hoc: on-current share of the nearly-free-electron band of [100] n = 3 (HW3 vs HW2) | `posthoc_nfe_share.py` | `fits/posthoc_2026-11-20/` |

The order above is the dependency order: `sweep` needs `hcalib`; `t1_cg.py` needs `rows/` and `dft.json`;
`t1_corr.py` needs the cache and the C_G results; `t1_pitch.py` and `t1_array_post.py` need `corr/`.

#### T1 band cache
`fits/t1/cache/` (46 `.npz` files, 143 MB) is not in the git repository. It is attached to the GitHub release
v1.0.0 as `t1_cache.tar.gz`: download it and run `tar -xzf t1_cache.tar.gz` in the repository root
(`sha256sum -c t1_cache.tar.gz.sha256` checks it). Or rebuild it:

```bash
PYTHONPATH=. python -P scripts/t1_cg.py cache     # one process, about 19 min, peak 2.04 GB (fits/t1/cg/README.md)
```

The build checks every ribbon against the parallel-plate I_ON stored in `fits/t1/rows/`. `t1_cg.py tob`,
`t1_corr.py run`, `t1_cg_array.py tob` and `t1_fig6.py` read the cache.

### 4. DFT-stage figures
`plot_a4b1.py` and `plot_a4b2.py` compare the calibration ribbons' Wannier bands with the DFT eigenvalues;
`render_structures.py` draws the bulk and sheet structures (OVITO environment). Fig. 1 is built in two steps:
`bash scripts/ovito_env.sh scripts/render_ribbons.py render` renders each panel with OVITO (set `OVITO_VENV` to the
OVITO environment; the wrapper runs it headless), then `PYTHONPATH=. python -P scripts/render_ribbons.py assemble`
lays out the panels in the analysis environment. `scripts/figstyle.py` holds the shared print style and the
journal's figure widths.

## Pre-registration

The fit recipes, pass criteria and device settings were written down and checksummed before the results they govern
were computed. The documents are deposited at the paths they had during the project:

| Document | Governs | Checksums |
|---|---|---|
| `fits/PREREG_sk_bulk_sheet.md` | joint bulk + sheet SK fit (gate a) | none |
| `fits/PREREG_sk_holdouts.md` | SK on the hold-outs (gates b–d) | `.md5` |
| `fits/PREREG_bb_holdouts.md` | building blocks on the hold-outs | `.md5`, `.code.md5` |
| `fits/PREREG_t1.md` | T1 width sweep, H-termination rule | `.md5`, `.code.md5` |
| `fits/PREREG_cg.md` | Laplace C_G | `.md5`, `.code.md5` |
| `fits/PREREG_corr.md` | Fröhlich B, tunnelling guard | `.md5`, `.code.md5` |
| `fits/PREREG_t1b.md` | pitch-normalised I_ON, error bands | `.md5`, `.code.md5` |
| `fits/PREREG_array.md` | periodic-array C_G | `.md5`, `.code.md5` |
| `dft/A3/README.md`, lines 1–84 | slab hold-out criteria | `dft/A3/README.prereg.md5` |

The DFT stage READMEs (`dft/*/README.md`) also carry each stage's registered criteria and amendments, together with the
run log. Each `.code.md5` lists the code the document froze; where a file was re-hashed later, the line comment says
why, and the last hash for a file is the one that applies. Those code files are deposited unchanged and match their
last hashes. The scripts that ran under a frozen document check its md5 and refuse to run if it has changed. Stored
results record the md5 of the document they ran under (`prereg_md5`).

The deposited copies of all eight `fits/PREREG_*.md` documents and of lines 1–84 of `dft/A3/README.md` differ from the
frozen originals in non-normative wording only: decision and approval notes are written as a neutral, dated record, a
remark on who made the recipe choices is removed or reworded, and pointers to project files that are not deposited are
removed or replaced by files in this repository. No criterion, number, date or rule changed. Each `.md5` file, and
`dft/A3/README.prereg.md5`, gives the deposited copy's hash first (the scripts check this one), followed by a note
with the original frozen hash. Results computed during the project record the original hash as `prereg_md5`; a rerun
records the deposited one. The rest of the stage READMEs was edited in the same way.

To check by hand from the repository root:

```bash
md5sum fits/PREREG_{sk_holdouts,bb_holdouts,t1,cg,corr,array}.md   # compare with the first hash in each .md5
(cd fits && md5sum PREREG_t1b.md)                                   # PREREG_t1b.md5 gives the path relative to fits/
sed -n 1,84p dft/A3/README.md | md5sum                              # compare with the first hash in dft/A3/README.prereg.md5
```

## Labels used in the code and notes

The docstrings, stage READMEs and pre-registration documents use the project's stage, gate and plan labels. The
project plan, its evidence record and the lab notebook are not deposited; this list gives what each label means here.

**Stages** (the DFT stages have their own READMEs under `dft/`):

| Label | Stage |
|---|---|
| T0 | laptop scaffold: the SK builder and fitter tested on digitised literature bands (`t0_sk_literature_fit.py`, `fits/T0_literature/`) |
| A0 | cutoff and pseudopotential tests (bulk ScN and the M-ScN sheet) |
| A1 | bulk rocksalt ScN, PBE and HSE06 (A1-1), then Wannier90 (A1-2) |
| A2 | the M-ScN sheet, PBE and HSE06 (A2-1), then Wannier90 (A2-2); `dft/A2/costtest/` holds the HSE cost tests |
| A3 | the three-plane slab hold-out: PAW PBE relaxation (A3-1), NC PBE twin and HSE06 (A3-2), Wannier90 (A3-3) |
| A4a | the four calibration ribbons: PBE relaxations (A4a-1); single points, bands and the magnetism, k and vacuum checks (A4a-2) |
| A4b | the calibration ribbons' edge Wannier recipe: NC PBE twin + Wannier90 (A4b-1), HSE06 + Wannier90 (A4b-2) |
| A5 | the two hold-out ribbons: PAW PBE relaxations (A5-1), magnetism and bands (A5-1b), NC PBE twin + Wannier90 (A5-2), HSE06 gap single points (A5-3) |
| T1 | laptop device stage: building-block ribbons at 1–6 nm, ToB metrics, and the C_G, Fröhlich and tunnelling corrections (`fits/t1/`) |
| T1b | the T1 pitch-normalised I_ON and E_g(W) error bands (`fits/PREREG_t1b.md`) |

Inside `dft/A2/costtest/` and in "cost test T1 / T2" elsewhere under `dft/`, T1 and T2 are the two HSE cost tests
(does smearing remove the empty-band EXX cost; how low the Wannier outer window can go, which set the ribbon rule
VBM + 11.05 eV in HSE06 and VBM + 10.23 eV in PBE), not the device stage.

**Gates**
- **Gate M** (`dft/A4a/README.md`): judged per calibration ribbon before any ribbon HSE06 run: magnetism, k
  convergence, PBE gap, edge states, vacuum.
- **Claims gates (a)–(e):** (a) the joint bulk + sheet SK fit (`fits/PREREG_sk_bulk_sheet.md`); (b) SK on the slab
  hold-out; (c) SK with its own edge blocks against the calibration ribbons' Wannier Hamiltonians (in-sample); (d) SK
  on the two ribbon hold-outs (`fits/PREREG_sk_holdouts.md`, `fits/holdouts/`); (e) the ToB code on the HfS2
  benchmark (`fits/gate_e/`).

**Wannier criteria and terms** (`dft/A4b/README.md`, `dft/README.md`)
- **P1–P6**, PBE twin, per ribbon and variant: P1 windows and band counts; P2 largest error in the device window
  ≤ 5 meV; P3 RMS in the frozen window ≤ 3 meV; P4 gap ≤ 3 meV; P5 band-edge masses ±3 %; P6 centres and spreads
  (descriptive).
- **H1–H4**, HSE06: H1 windows; H2 ≤ 1 meV against the HSE grid eigenvalues in the frozen window; H3 a stop flag if a
  non-model state lies in the device window; H4 (added 2026-11-14) every pair of variants within 5 meV over the
  device window.
- **Device window:** [VBM − 0.5 eV, CBM + 0.5 eV]. **Non-model state:** an empty state with s_out ≥ 0.25, i.e. weight
  largely outside the projected basis.
- **W1, W2, W3, HW1–HW3; sets B and S:** the Wannier variants and projection sets (`dft/README.md`).
- **Route O / route N:** HSE06 on the irreducible grid unfolded with `open_grid.x`, or on the full grid without
  symmetry (`dft/A1/README.md`).
- **PBE twin:** the norm-conserving PBE counterpart of an HSE06 calculation, Wannier recipe or SK fit; the hold-out
  transfer tests compare PBE twins.

**Plan labels**
- **§3, §4, §10:** the plan's design items (D labels), its stages and claims gates, and its risk register (LR labels).
- **D3:** the SK model specification: Sc d + N p basis, shells up to the third-neighbour Sc–Sc, Harrison scaling,
  surface on-site blocks, our own NC-HSE06 fit targets, the projected Wannier H(R) as the initial guess, and the
  thickness axis.
- **D4:** the ribbon construction: orbital-resolved edge blocks read from the calibration ribbons' Wannier H(R), the
  calibration edge relaxation copied onto wider ribbons, and the E_g(W) error bands.
- **D6:** the device corrections: geometric (Laplace) C_G, Fröhlich ballisticity, tunnelling guard, pitch-normalised
  I_ON.
- **D7:** the transfer tests on the hold-outs and their pass criteria (gap and valley offsets ±50 meV, masses ±10 %,
  I_ON ±15 %, the same IRDS verdict; SK-HSE06 gaps ±80 meV).
- **LR1:** the fallback if gate (a) fails: DFT-Wannier device Hamiltonians for the ribbons, or the SK errors carried
  into the error bars. Its option A (HSE06 + Wannier90 device Hamiltonians for the five DFT ribbons, `dft/A4a/cost_a4.py`)
  was later reversed for the two hold-outs (`dft/A5/README.md`); option B adds two wider ribbons ([110] n = 11,
  [100] n = 8).
- **LR2:** the fallback if the hold-out gates fail: no refitting, the failure reported as a transferability limit, and
  building-block Wannier ribbons as the device Hamiltonians.
- **LR3:** an edge state sets the gap: characterise it, make the H-terminated ribbon the device case and report both
  edge chemistries.
- **LR4:** a valley crossover outside the DFT bracket is reported as a model prediction (triggered by an M̄ valley
  offset > 0.6 eV).
- **LR5:** edge magnetism (ΔE_mag > 5 meV per edge atom): a spin-resolved edge term, or H-terminated edges only.
- **LR6:** no Wannier variant passes: recount the bands and narrow the window.
- **LR10:** the HfS2 benchmark misses ±10 %: flagged with its suspected cause; no convention changes.
- **LK7, LK11, LM7, LM8, K3, F1, F2, M16:** items of the planning evidence (its kill panel, referee items and cost
  tables), cited where a requirement came from. In context: LK7, the slab hold-out along the thickness axis; LK11, the
  ballisticity bound and the pitch-normalised comparison; LM7, the ≈ 0.13 eV semicore gap shift between literature
  HSE06 sets; LM8, the Fröhlich mean free path; K3, the geometric C_G by 2D Laplace; F1, the `open_grid.x` route
  saving; F2, the cost estimate and the cost law used in `dft/A4a/cost_a4.py`; M16 ("heavy-plan M16"), the
  tunnelling guard's λ_s range and R > 1.5 flag. The "heavy plan" is an earlier, DFT-heavier version of the plan.
- **The notebook, `00_Master_Notebook.md`:** the project's lab notebook. Its "Pre-registered device conventions"
  (signed off 2026-11-11) are restated in the `scnfet/tob.py` docstring and in the device settings of the
  pre-registration documents.

## HfS2 benchmark data

`gate_e_hfs2_benchmark.py` reproduces the ToB benchmark against Matić and Poljak, Electronics 13, 1048 (2024). Its
input Hamiltonians, QUDEN_MATERIALS (Poljak and Matić; CC BY-NC-SA 4.0), are not redistributed here. Download them
from the authors' distribution and place the `HfS2XNR` set at
`benchmark/quden/materials/MATERIALS_QUDEN/HfS2XNR/<N>HfS2XNR_hr.dat`. The stored results are in `fits/gate_e/`.

## What is not included

The git repository holds everything except `fits/t1/cache/`, which is a release asset ("T1 band cache" above).
Neither holds wavefunctions, charge densities, potential cube files, Wannier90 overlap and projection files (`.amn`,
`.mmn`, `.chk`, `.eig`), `_wsvec.dat` files, or the raw QE outputs that no script reads. The Wannier Hamiltonians included
are those the scripts read (see `dft/README.md`); pseudopotentials are the published SSSP, SG15 and PseudoDojo files
named in the inputs.

## Licence

Code (`scnfet/`, `scripts/`, the `*.py` and `*.sh` files under `dft/`): MIT, see `LICENSE`.
Data (everything else, including `dft/` inputs and outputs and `fits/`): CC BY 4.0, see `LICENSE-DATA`.

## Citation

See `CITATION.cff`. Please cite the paper as well once it is published.
