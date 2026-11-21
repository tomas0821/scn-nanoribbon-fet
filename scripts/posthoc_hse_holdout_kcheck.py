"""Post hoc, descriptive, not pre-registered (2026-11-19; Supplementary Section S1): HSE06 building-block bands of the two bare hold-outs at the
A5-3 HSE06 SCF k-points, against the DFT HSE06 eigenvalues there, over the device window. Read-only."""
import importlib.util, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from posthoc_qe_parse import last_block

ROOT = Path('.').resolve()
spec = importlib.util.spec_from_file_location('gbb', ROOT / 'scripts' / 'gate_bb_holdouts.py')
GB = importlib.util.module_from_spec(spec); argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(GB); sys.argv = argv
H2M = 3.80998212
for kh, kc in GB.HOLD.items():
    st = GB.MBH.geometry(kh); nv = 3 * st.species.count('N')
    hr, wf, used, info = GB.build(GB.MBH, kh, kc, 'HSE06')
    ks, nocc, a = last_block(str(ROOT / f'dft/A5/hse/outputs/hse_{kh}.out'))
    kx = np.array([abs(k) for k, _ in ks]); Ed = np.array([e for _, e in ks])
    kf = np.c_[kx, np.zeros((len(kx), 2))]
    Eb = GB.G.w_eigs(hr, kf)
    # DFT band j (0-based) ↔ BB band j - (nocc - nv)
    off = nocc - nv
    vd, cd = Ed[:, nocc - 1].max(), Ed[:, nocc].min()
    vb, cb = Eb[:, nv - 1].max(), Eb[:, nv].min()
    print(f'\n{kh}: BB HSE06 {info["variant"]}, shift {info["sheet_shift_eV"]:+.3f} eV; gap BB {cb-vb:.4f} vs DFT {cd-vd:.4f} '
          f'({1e3*((cb-vb)-(cd-vd)):+.1f} meV); a = {a:.3f} Å, {len(kx)} k')
    worst_v = worst_c = worst_raw = 0.0
    for ik in range(len(kx)):
        lines = []
        for j in range(nocc - 6, nocc + 4):
            ed = Ed[ik, j]; eb = Eb[ik, j - off]
            if j < nocc and ed >= vd - 0.5:
                d = 1e3 * ((eb - vb) - (ed - vd)); worst_v = max(worst_v, abs(d)); lines.append(f'VB{j-nocc+1:+d} {ed-vd:+.3f}:{d:+.1f}')
            if j >= nocc and ed <= cd + 0.5:
                d = 1e3 * ((eb - cb) - (ed - cd)); dr = 1e3 * ((eb - vb) - (ed - vd))
                worst_c = max(worst_c, abs(d)); worst_raw = max(worst_raw, abs(dr))
                lines.append(f'CB{j-nocc:+d} {ed-cd:+.3f}:{d:+.1f}')
        if lines:
            print(f'  k {kx[ik]:.4f}: ' + '  '.join(lines) + '   (E_DFT rel. edge, eV : BB−DFT, meV, edge-aligned)')
    print(f'  max |dev| in window, each side edge-aligned: VB {worst_v:.1f} meV, CB {worst_c:.1f} meV; '
          f'CB with VBM alignment (incl. gap error) {worst_raw:.1f} meV')
    k1 = kx[1] * 2 * np.pi / a
    for name, jd in (('CB', nocc), ('VB', nocc - 1)):
        md = H2M * k1 ** 2 / (Ed[1, jd] - Ed[0, jd]); mb = H2M * k1 ** 2 / (Eb[1, jd - off] - Eb[0, jd - off])
        print(f'  {name} chord mass Γ→k1 ({abs(Ed[1, jd]-Ed[0, jd])*1e3:.0f} meV): DFT {md:+.3f}, BB {mb:+.3f} ({100*(mb/md-1):+.1f} %)')
