"""Post hoc (2026-11-19; Supplementary Section S1): the same Γ→k1 chord comparison in PBE (BB PBE vs the PBE-Wannier W2 hold-out), to calibrate the HSE06
chord numbers against the registered PBE test. Read-only."""
import importlib.util, sys, json
from pathlib import Path
import numpy as np
ROOT = Path('.').resolve()
spec = importlib.util.spec_from_file_location('gbb', ROOT / 'scripts' / 'gate_bb_holdouts.py')
GB = importlib.util.module_from_spec(spec); argv = sys.argv; sys.argv = ['x']; spec.loader.exec_module(GB); sys.argv = argv
W = GB.W
H2M = 3.80998212
for kh, kc, dk in (('r110_8', 'r110_4', 0.0625), ('r100_6', 'r100_3', 1 / 12)):
    st = GB.MBH.geometry(kh); nv = 3 * st.species.count('N'); a = st.cell[0, 0]
    hr_w = W.read_hr(GB.M5.TWIN / 'w90' / kh / 'W2' / f'{kh}_W2_hr.dat')
    hr_b, wf_b, _, _ = GB.build(GB.MBH, kh, kc, 'PBE')
    kf = np.array([[0, 0, 0], [dk, 0, 0], [2 * dk, 0, 0], [3 * dk, 0, 0]])
    Ew, Eb = GB.G.w_eigs(hr_w, kf), GB.G.w_eigs(hr_b, kf)
    k1 = dk * 2 * np.pi / a
    for name, j in (('CB', nv), ('VB', nv - 1)):
        mw = H2M * k1 ** 2 / (Ew[1, j] - Ew[0, j]); mb = H2M * k1 ** 2 / (Eb[1, j] - Eb[0, j])
        print(f'{kh} PBE {name} chord Γ→k1 ({abs(Ew[1, j]-Ew[0, j])*1e3:.0f} meV): Wannier {mw:+.3f}, BB {mb:+.3f} ({100*(mb/mw-1):+.1f} %);'
              f'  ΔE(k1..k3) BB−W edge-aligned (meV): ' + ' '.join(f'{1e3*((Eb[i,j]-Eb[0,j])-(Ew[i,j]-Ew[0,j])):+.1f}' for i in (1, 2, 3)))
    gb = json.loads((ROOT / 'fits/bb/gate_bb.json').read_text())['ribbons'][kh]
    print('   registered PBE masses: W2', {k: round(gb['wannier'][k], 3) for k in ('m_VB', 'm_CB')},
          'BB', {k: round(gb['BB PBE']['metrics'][k], 3) for k in ('m_VB', 'm_CB')},
          'BB HSE06', {k: round(gb['BB HSE06']['metrics'][k], 3) for k in ('m_VB', 'm_CB')})
