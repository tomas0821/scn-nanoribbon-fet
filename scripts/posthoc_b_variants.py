"""Post hoc sensitivity (2026-11-19, not pre-registered; Supplementary Section S4, Table S2): recompute B at eps12 = 1, 4, 5 and with lambda0/(lambda0+2L); recount IRDS passes per W, EOT 0.6.
Read-only on project data. Uses stored ballistic I_ON (fits/t1/corr/*.json, fits/t1/cg_array/array.json)."""
import json, sys
from pathlib import Path
from scnfet import frohlich as FR

ROOT = Path('.')
IRDS, L = 1.979, 15.0
FAMS = ('100', '110', '110h')
EPS = (1.0, 4.0, 5.0)
arr = json.loads((ROOT / 'fits/t1/cg_array/array.json').read_text())['families']
rows = {}
for fam in FAMS:
    a = {r['n']: r for r in arr[fam]}
    for p in sorted((ROOT / 'fits/t1/corr').glob(f'{fam}_n*.json'), key=lambda p: int(p.stem.split('_n')[1])):
        d = json.loads(p.read_text()); n = d['n']
        for c in ('n', 'p'):
            r = d['carriers'][c]; m = r['m']
            lam = {}
            for e in EPS:
                lam[e] = FR.Frohlich(m, eps12=e).lam0('2d')
            lam3 = r['frohlich']['3d']['lam0_nm']
            assert abs(lam[1.0] / r['frohlich']['2d']['lam0_nm'] - 1) < 1e-6, (fam, n, c)
            I = {g: r['I_on'][f'{g} EOT 0.6']['ballistic'] for g in ('parallel_plate', 'single', 'gaa')}
            I['s2'] = a[n]['s2'][c]['per_W']['ballistic']; I['s7.7'] = a[n]['sGAA'][c]['per_W']['ballistic']
            rows[(fam, n, c)] = dict(W=d['W_nm'], m=m, lam=lam, lam3=lam3, I=I)

def B(l, LL=L):
    return l / (l + LL)

# B ranges
print('== B ranges (2D eps12 = 1/4/5; 3D) ==')
for fam in FAMS:
    for c in ('n', 'p'):
        rr = [v for k, v in rows.items() if k[0] == fam and k[2] == c and v['m'] < 5]
        s = ' '.join(f"e{e:g}: {min(B(v['lam'][e]) for v in rr):.3f}-{max(B(v['lam'][e]) for v in rr):.3f}" for e in EPS)
        s3 = f"3D: {min(B(v['lam3']) for v in rr):.3f}-{max(B(v['lam3']) for v in rr):.3f}"
        s2L = f"2D e1 2L: {min(B(v['lam'][1.0], 2*L) for v in rr):.3f}-{max(B(v['lam'][1.0], 2*L) for v in rr):.3f}"
        print(f'{fam:5s} {c}: m {min(v["m"] for v in rr):.2f}-{max(v["m"] for v in rr):.2f} | {s} | {s3} | {s2L}')

GE = ('parallel_plate', 's2', 's7.7', 'single', 'gaa')
def count(fun):
    out = {}
    for fam in FAMS:
        for c in ('n', 'p'):
            ks = sorted([k for k in rows if k[0] == fam and k[2] == c], key=lambda k: k[1])
            out[(fam, c)] = {g: (sum(fun(rows[k], g) >= IRDS for k in ks), len(ks),
                                 [round(rows[k]['W'], 2) for k in ks if fun(rows[k], g) >= IRDS]) for g in GE}
    return out

cases = {
    'ballistic': lambda v, g: v['I'][g],
    '2D e1 (registered)': lambda v, g: v['I'][g] * B(v['lam'][1.0]),
    '3D': lambda v, g: v['I'][g] * B(v['lam3']),
    '2D e4': lambda v, g: v['I'][g] * B(v['lam'][4.0]),
    '2D e5': lambda v, g: v['I'][g] * B(v['lam'][5.0]),
    '2D e1, 2L': lambda v, g: v['I'][g] * B(v['lam'][1.0], 2 * L),
    '3D, 2L': lambda v, g: v['I'][g] * B(v['lam3'], 2 * L),
    '2D e5, 2L': lambda v, g: v['I'][g] * B(v['lam'][5.0], 2 * L),
}
for name, fun in cases.items():
    print(f'\n== pass counts per W, EOT 0.6: {name} ==  (pp / s2 / s7.7 / single / gaa)')
    for (fam, c), d in count(fun).items():
        print(f'{fam:5s} {c}: ' + ' / '.join(f'{d[g][0]}/{d[g][1]}' for g in GE))
        if name in ('2D e5', '2D e1, 2L', '2D e4'):
            for g in GE:
                if 0 < d[g][0] < d[g][1]:
                    print(f'        {g}: passing W = {d[g][2]}')
# closest-to-line listing for the registered and 3D cases
print('\n== margins near the line (|I/IRDS - 1| < 2%) ==')
for name in ('2D e1 (registered)', '3D'):
    fun = cases[name]
    for k, v in sorted(rows.items()):
        for g in GE:
            x = fun(v, g) / IRDS - 1
            if abs(x) < 0.02:
                print(f'{name:20s} {k[0]:5s} n={k[1]:2d} {k[2]} W {v["W"]:.2f} {g:15s} I {fun(v, g):.4f}  {100*x:+.2f}%')
