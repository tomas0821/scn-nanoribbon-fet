"""Parse the last eigenvalue block of a QE pw.x output; count bands in the device window per k."""
import re, sys
import numpy as np

def last_block(path):
    txt = open(path).read()
    i = txt.rfind('End of self-consistent calculation')
    blk = txt[i:]
    ks = []
    for m in re.finditer(r'k =\s*([-\d.]+)\s*([-\d.]+)\s*([-\d.]+).*?bands \(ev\):\s*\n(.*?)(?:occupation numbers)', blk, re.S):
        e = np.array([float(x) for x in re.findall(r'-?\d+\.\d+', m.group(4))])
        ks.append((float(m.group(1)), e))
    nel = float(re.search(r'number of electrons\s*=\s*([\d.]+)', txt).group(1))
    nalat = float(re.search(r'lattice parameter \(alat\)\s*=\s*([\d.]+)', txt).group(1)) * 0.529177
    return ks, int(round(nel / 2)), nalat

for p in sys.argv[1:]:
    ks, nocc, a = last_block(p)
    E = np.array([e for _, e in ks])
    vbm, cbm = E[:, nocc - 1].max(), E[:, nocc].min()
    print(f'{p.split("/")[-1]}: {len(ks)} k (Γ→X), a = {a:.3f} Å, Δk = {abs(ks[1][0]-ks[0][0]):.4f} 2π/a '
          f'= {abs(ks[1][0]-ks[0][0])*2*np.pi/a:.4f} 1/Å; nbnd {E.shape[1]}, nocc {nocc}; '
          f'VBM {vbm:.4f} CBM {cbm:.4f} gap {cbm-vbm:.4f} eV')
    for kx, e in ks:
        nv = int(((e[:nocc] >= vbm - 0.5)).sum()); nc = int(((e[nocc:] <= cbm + 0.5)).sum())
        print(f'   k {kx:+.4f}: VB in [VBM-0.5,VBM] {nv:2d}  CB in [CBM,CBM+0.5] {nc:2d}   '
              f'VB top {e[nocc-1]-vbm:+.4f}  CB bottom {e[nocc]-cbm:+.4f}  CB2 {e[nocc+1]-cbm:+.4f}  VB2 {e[nocc-2]-vbm:+.4f}')
    # finite-difference parabolic masses from the first two k-points of the edge bands
    H2M = 3.80998212
    k1 = abs(ks[1][0]-ks[0][0])*2*np.pi/a
    for name, j in (('CB', nocc), ('VB', nocc - 1)):
        dE = E[1, j] - E[0, j]
        print(f'   {name} at Γ→k1: ΔE {1e3*dE:+.1f} meV, m* (parabola through Γ, k1) {H2M*k1**2/dE:+.3f}')
