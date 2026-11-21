"""Post hoc (2026-11-19; the dagger marks of Table 3): Table 3 cells — widths within ±4 % (and ±8 %) of the IRDS line, from the stored values. Read-only."""
import json
from pathlib import Path
IRDS = 1.979
arr = json.loads(Path('fits/t1/cg_array/array.json').read_text())['families']
cells = []
for fam in ('100', '110', '110h'):
    a = {r['n']: r for r in arr[fam]}
    pts = []
    for p in sorted(Path('fits/t1/corr').glob(f'{fam}_n*.json'), key=lambda p: int(p.stem.split('_n')[1])):
        d = json.loads(p.read_text()); n = d['n']
        for c in ('n', 'p'):
            I = d['carriers'][c]['I_on']
            v = {('pp', 'ball'): I['parallel_plate EOT 0.6']['ballistic']}
            for g, key in (('pp', 'parallel_plate'), ('single', 'single'), ('gaa', 'gaa')):
                for b in ('B2d', 'B3d'):
                    v[(g, b)] = I[f'{key} EOT 0.6'][b]
            for g, key in (('s2', 's2'), ('s7.7', 'sGAA')):
                for b in ('B2d', 'B3d'):
                    v[(g, b)] = a[n][key][c]['per_W'][b]
            pts.append((c, n, d['W_nm'], v))
    for c in ('n', 'p'):
        for col in [('pp', 'ball')] + [(g, b) for g in ('pp', 's2', 's7.7', 'single', 'gaa') for b in ('B2d', 'B3d')]:
            xs = [(n, W, v[col] / IRDS - 1) for cc, n, W, v in pts if cc == c]
            near4 = [(round(W, 2), round(100 * x, 2)) for n, W, x in xs if abs(x) < 0.04]
            near8 = sum(abs(x) < 0.08 for n, W, x in xs)
            npass = sum(x >= 0 for n, W, x in xs)
            if near4:
                print(f'{fam:5s} {c} {col[0]:6s} {col[1]:4s}: pass {npass:2d}/{len(xs)}; within ±4 %: {len(near4)} {near4}; within ±8 %: {near8}')
