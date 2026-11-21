"""Fig. 1: structures. Bulk ScN, the M-ScN sheet and the DFT ribbons at their relaxed PBE geometries, with W marked.

Two stages (OVITO and matplotlib live in different environments):

    bash scripts/ovito_env.sh scripts/render_ribbons.py render       # ~/.venvs/ovito (OVITO 3.17, Tachyon)
    PYTHONPATH=. ~/.venvs/scnfet/bin/python -P scripts/render_ribbons.py assemble

`render` writes the panel images and their geometry to figures/fig1_renders/ (PNG + meta.json); `assemble` lays
them out as figures/Fig1_ribbons.png at the sn-jnl text width (5.147 in; scripts/figstyle.py).

Panels (3 × 3):
  (a) bulk rocksalt ScN, a = 4.501 Å (the A1 geometry, experimental a);
  (b) the M-ScN sheet (the A2 geometry: A0 PBE sheet scaled by a_exp / a_PBE; a_p 3.109 Å, planes 2.17 Å apart);
  (c) bare [110] n = 4 seen along its axis (both planes; relaxed);
  (d–f) [100]: n = 3 bare, n = 3 H (calibration ribbons), n = 6 (hold-out);
  (g–i) [110]: n = 4 bare, n = 4 H (calibration ribbons), n = 8 (hold-out).
Ribbon geometries are the final BFGS coordinates of the PBE relaxations (dft/A4a/outputs/relax_*.out,
dft/A5/outputs/relax_*.out). Top views (d–i) share one scale and show the upper plane; H atoms of the lower plane
are drawn faded, with their bonds, because at some edges they are the outermost nuclei. The lower plane's host atoms
lie beneath the upper plane's (other species) within 0.05 Å in y. W is the edge-to-edge width between the outermost
nuclei of both planes, H included (Methods), i.e. scnfet.structures.width of the relaxed structure; the printed
values are the relaxed ones of Table 1.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "Fig1_ribbons.png"
RDIR = ROOT / "figures" / "fig1_renders"
BOHR = 0.529177210903
COL = {"Sc": "#c9b48c", "N": "#4a64d6", "H": "#e9e9e9"}
RAD = {"Sc": 0.56, "N": 0.40, "H": 0.30}
BOND_CUT = {("N", "Sc"): 2.6, ("H", "Sc"): 2.1, ("H", "N"): 1.3}
LOWER_ALPHA = 0.62                       # transparency of the lower-plane H and their bonds (top views)
RIBBONS = {  # key: (relax output, direction, label)
    "r100_3": ("dft/A4a/outputs/relax_r100_3.out", "100", "[100] $n$ = 3"),
    "r100_3h": ("dft/A4a/outputs/relax_r100_3h.out", "100", "[100] $n$ = 3 H"),
    "r100_6": ("dft/A5/outputs/relax_r100_6.out", "100", "[100] $n$ = 6, hold-out"),
    "r110_4": ("dft/A4a/outputs/relax_r110_4.out", "110", "[110] $n$ = 4"),
    "r110_4h": ("dft/A4a/outputs/relax_r110_4h.out", "110", "[110] $n$ = 4 H"),
    "r110_8": ("dft/A5/outputs/relax_r110_8.out", "110", "[110] $n$ = 8, hold-out"),
}
TABLE1_W_NM = {"r100_3": 1.09, "r100_3h": 1.32, "r110_4": 1.07, "r110_4h": 1.39, "r110_8": 2.31, "r100_6": 2.41}
# top views: one scale; visible window (Å)
WIN_X, WIN_Y = 21.0, 26.4
PX_PER_A = 46                            # render resolution (≈ 2× the printed 300 dpi at the figure's scale)
# A1 bulk and A2 sheet geometry
A_BULK = 4.501
A2_A_PBE, A2_C, A2_SCALE = 3.117899, 22.5, 4.501 / 4.5134
A2_Z = {"Sc1": 0.4516957094, "N1": 0.4516288717, "N2": 0.5483711283, "Sc2": 0.5483042906}


# ================================================================================================= geometry
def relaxed(path):
    """(species, positions Å, cell rows Å) from the final coordinates of a QE relax output."""
    txt = (ROOT / path).read_text()
    assert "bfgs converged" in txt and "End final coordinates" in txt, f"{path}: not converged"
    blk = txt[txt.rindex("Begin final coordinates"):txt.rindex("End final coordinates")]
    blk = blk.split("ATOMIC_POSITIONS")[1]
    assert blk.lstrip().startswith("(angstrom)"), "expected angstrom positions"
    rows = re.findall(r"^\s*([A-Za-z]+)\d*\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", blk, re.M)
    alat = float(re.search(r"lattice parameter \(alat\)\s+=\s+([\d.]+)", txt).group(1)) * BOHR
    cell = np.array([[float(v) for v in re.search(rf"a\({i}\) = \(\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", txt).groups()]
                     for i in (1, 2, 3)]) * alat
    return [r[0] for r in rows], np.array([[float(v) for v in r[1:]] for r in rows]), cell


def width(pos):
    return float(np.ptp(pos[:, 1]))


def planes(sp, pos):
    """Boolean mask of the upper plane (z above the host atoms' mean; each H goes with its nearest host atom)."""
    host = np.array([s != "H" for s in sp])
    zmid = pos[host, 2].mean()
    up = pos[:, 2] > zmid
    for i in np.where(~host)[0]:
        d = np.linalg.norm(pos[host] - pos[i], axis=1)
        up[i] = pos[host][np.argmin(d), 2] > zmid
    return up


def bonds(sp, pos):
    out = []
    for i in range(len(sp)):
        for j in range(i + 1, len(sp)):
            c = BOND_CUT.get(tuple(sorted((sp[i], sp[j]))))
            if c and np.linalg.norm(pos[i] - pos[j]) < c:
                out.append((i, j))
    return np.array(out, int).reshape(-1, 2)


def repeat_x(sp, pos, L, x0, x1):
    """Periodic images along x covering [x0, x1]."""
    m0, m1 = int(np.floor((x0 - pos[:, 0].max()) / L)), int(np.ceil((x1 - pos[:, 0].min()) / L))
    P = np.vstack([pos + [k * L, 0, 0] for k in range(m0, m1 + 1)])
    return list(sp) * (m1 - m0 + 1), P


def bulk_cube(a=A_BULK):
    pos, sym = [], []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                pos.append(np.array([i, j, k]) * a / 2)
                sym.append("Sc" if (i + j + k) % 2 == 0 else "N")
    return sym, np.array(pos)


def sheet_patch(n=5):
    """A2 sheet on the conventional (Sc–N bond) grid, both planes: lower plane Sc1/N1, upper N2/Sc2."""
    a_c = A2_A_PBE * A2_SCALE * np.sqrt(2)
    z = {k: (v - 0.5) * A2_C * A2_SCALE for k, v in A2_Z.items()}
    pos, sym = [], []
    for p in range(2):
        for i in range(n + 1):
            for j in range(n + 1):
                s = "Sc" if (i + j + p) % 2 == 0 else "N"
                zz = (z["Sc1"] if s == "Sc" else z["N1"]) if p == 0 else (z["Sc2"] if s == "Sc" else z["N2"])
                pos.append([i * a_c / 2, j * a_c / 2, zz])
                sym.append(s)
    return sym, np.array(pos), a_c


# ================================================================================================= render (OVITO)
def _render(sym, pos, transp, bnd, btransp, direction, up, center, win, fname, cell=None):
    from ovito.data import DataCollection, ParticleType
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import Viewport, TachyonRenderer
    d = DataCollection()
    if cell is None:
        c = d.create_cell(np.c_[np.diag([300.0] * 3), np.full(3, -150.0)], pbc=(False, False, False))
        c.vis.enabled = False
    else:
        c = d.create_cell(np.c_[np.asarray(cell, float).T, np.zeros(3)], pbc=(False, False, False))
        c.vis.line_width = 0.07
        c.vis.rendering_color = (0.33, 0.35, 0.37)
    p = d.create_particles(count=len(sym))
    p.create_property("Position", data=pos)
    tp = p.create_property("Particle Type")
    names = sorted(set(sym), key=["Sc", "N", "H"].index)
    for k, el in enumerate(names):
        rgb = tuple(int(COL[el][i:i + 2], 16) / 255 for i in (1, 3, 5))
        tp.types.append(ParticleType(id=k + 1, name=el, color=rgb, radius=RAD[el]))
    tp[...] = [names.index(s) + 1 for s in sym]
    p.create_property("Transparency", data=np.asarray(transp, float))
    if len(bnd):
        b = p.create_bonds(count=len(bnd))
        b.create_property("Topology", data=bnd)
        b.create_property("Transparency", data=np.asarray(btransp, float))
        b.create_property("Color", data=np.tile([0.60, 0.61, 0.63], (len(bnd), 1)))
        b.vis.width = 0.15
    pipe = Pipeline(source=StaticSource(data=d))
    pipe.add_to_scene()
    dvec = np.asarray(direction, float) / np.linalg.norm(direction)
    vp = Viewport(type=Viewport.Type.Ortho, camera_dir=tuple(dvec), camera_pos=tuple(np.asarray(center) - 150 * dvec),
                  fov=win[1] / 2)                         # ortho fov = half the vertical extent (checked by test render)
    vp.camera_up = tuple(up)
    size = (int(round(win[0] * PX_PER_A)), int(round(win[1] * PX_PER_A)))
    rend = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_brightness=0.85, ambient_occlusion_samples=12,
                           antialiasing_samples=9, shadows=False, direct_light_intensity=0.9)
    vp.render_image(size=size, filename=str(RDIR / fname), background=(1, 1, 1), renderer=rend)
    pipe.remove_from_scene()
    return size


def render():
    RDIR.mkdir(parents=True, exist_ok=True)
    meta = {"px_per_A": PX_PER_A, "win_top": [WIN_X, WIN_Y], "panels": {}}
    # ---- ribbons, top views (upper plane + faded lower-plane H)
    for key, (path, dr, label) in RIBBONS.items():
        sp, pos, cell = relaxed(path)
        W = width(pos)
        assert abs(W / 10 - TABLE1_W_NM[key]) < 0.005, (key, W)
        up = planes(sp, pos)
        L = cell[0, 0]
        yc = (pos[:, 1].min() + pos[:, 1].max()) / 2
        xc = pos[:, 0].mean()
        S, P = repeat_x(sp, pos, L, xc - WIN_X / 2 - 3, xc + WIN_X / 2 + 3)
        UP = np.tile(up, len(S) // len(sp))
        # lower-plane parents of lower-plane H: kept invisible so the faded bond has an anchor
        bnd_all = bonds(S, P)
        low_h = (~UP) & (np.array(S) == "H")
        anchor = np.zeros(len(S), bool)
        for i, j in bnd_all:
            if low_h[i] and not UP[j]:
                anchor[j] = True
            if low_h[j] and not UP[i]:
                anchor[i] = True
        keep = UP | low_h | anchor
        idx = np.where(keep)[0]
        remap = -np.ones(len(S), int)
        remap[idx] = np.arange(len(idx))
        bnd = np.array([(remap[i], remap[j]) for i, j in bnd_all if keep[i] and keep[j]], int).reshape(-1, 2)
        transp = np.where(UP[idx], 0.0, np.where(anchor[idx], 1.0, LOWER_ALPHA))
        btransp = np.array([0.0 if UP[idx[i]] and UP[idx[j]] else LOWER_ALPHA for i, j in bnd])
        zc = pos[:, 2].max()
        size = _render([S[k] for k in idx], P[idx], transp, bnd, btransp, (0, 0, -1), (0, 1, 0), (xc, yc, zc),
                       (WIN_X, WIN_Y), f"{key}.png")
        dz_host = float(np.abs(pos[(~up) & (np.array(sp) != "H"), 1].min() - pos[up & (np.array(sp) != "H"), 1].min()))
        meta["panels"][key] = dict(label=label, direction=dr, W_A=W, y_lo=float(pos[:, 1].min() - yc),
                                   y_hi=float(pos[:, 1].max() - yc), size=size, win=[WIN_X, WIN_Y],
                                   upper_plane_W_A=width(pos[up]), host_plane_edge_offset_A=dz_host,
                                   n_atoms=len(sp), source=path)
        print(f"{key}: W {W / 10:.3f} nm (upper plane {width(pos[up]) / 10:.3f}); render {size}", flush=True)
    # ---- (c) bare [110] n = 4 along its axis: both planes, three cells (they project onto each other)
    sp, pos, cell = relaxed(RIBBONS["r110_4"][0])
    S, P = repeat_x(sp, pos, cell[0, 0], pos[:, 0].min() - 0.1, pos[:, 0].max() + cell[0, 0])
    bnd = bonds(S, P)
    yc, zc = (P[:, 1].min() + P[:, 1].max()) / 2, P[:, 2].mean()
    win_c = (12.6, 5.6)
    size = _render(S, P, np.zeros(len(S)), bnd, np.zeros(len(bnd)), (-1, 0, 0), (0, 0, 1), (P[:, 0].mean(), yc, zc),
                   win_c, "r110_4_axis.png")
    zs = sorted({round(z, 1) for z in pos[:, 2]})
    meta["panels"]["r110_4_axis"] = dict(size=size, win=list(win_c), W_A=width(pos), y_lo=float(pos[:, 1].min() - yc),
                                         y_hi=float(pos[:, 1].max() - yc),
                                         z_planes=[float(np.mean(pos[pos[:, 2] < zc, 2]) - zc),
                                                   float(np.mean(pos[pos[:, 2] > zc, 2]) - zc)])
    print("r110_4_axis:", size, zs, flush=True)
    # ---- (a) bulk cube, oblique
    sym, pos = bulk_cube()
    size = _render(sym, pos, np.zeros(len(sym)), bonds(sym, pos), np.zeros(len(bonds(sym, pos))), (-1.0, -1.45, -0.75),
                   (0, 0, 1), np.full(3, A_BULK / 2), (9.4, 9.4), "bulk.png", cell=np.eye(3) * A_BULK)
    meta["panels"]["bulk"] = dict(size=size, a_A=A_BULK)
    print("bulk:", size, flush=True)
    # ---- (b) sheet, oblique
    sym, pos, a_c = sheet_patch(5)
    bd = bonds(sym, pos)
    c = [pos[:, 0].mean(), pos[:, 1].mean(), pos[:, 2].mean()]
    size = _render(sym, pos, np.zeros(len(sym)), bd, np.zeros(len(bd)), (-1.0, -1.55, -0.55), (0, 0, 1), c,
                   (17.6, 14.0), "sheet.png")
    meta["panels"]["sheet"] = dict(size=size, a_p_A=a_c / np.sqrt(2),
                                   plane_gap_A=float(np.ptp(pos[:, 2])))
    print("sheet:", size, f"a_p {a_c / np.sqrt(2):.4f}", flush=True)
    (RDIR / "meta.json").write_text(json.dumps(meta, indent=1))
    print("wrote", RDIR)


# ================================================================================================= assemble (matplotlib)
def assemble():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import figstyle as FS
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from PIL import Image

    FS.use(7.5)
    plt.rcParams["figure.constrained_layout.use"] = False
    meta = json.loads((RDIR / "meta.json").read_text())
    FW = FS.width(1.0)
    k = 0.0640                                   # in per Å, top views (one scale)
    col_w = FW / 3
    strip = 4.4                                  # Å right of each top view for the W arrow and label
    top_w, top_h = (WIN_X + strip) * k, WIN_Y * k
    title_h, gap_h, row1_h, leg_h = 0.17, 0.05, 1.30, 0.20
    FH = row1_h + title_h + 2 * (top_h + title_h + gap_h) + leg_h + 0.03
    fig = plt.figure(figsize=(FW, FH))

    def axes_in(x, y, w, h):                     # inches from the top-left corner
        return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH])

    def title(x, y, s, lab):                       # inches from the top-left corner (baseline just above y)
        fig.text(x / FW, 1 - y / FH, rf"$\mathbf{{({s})}}$ {lab}", ha="left", va="bottom", fontsize=8.0,
                 color=FS.INK)

    def blank(ax):
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    # ---- row 1: bulk, sheet, axial view
    y0 = title_h
    for j, (key, s, lab) in enumerate((("bulk", "a", "bulk ScN (rocksalt)"), ("sheet", "b", "M-ScN sheet"),
                                       ("r110_4_axis", "c", "[110] $n$ = 4, axial view"))):
        img = np.asarray(Image.open(RDIR / f"{key}.png").convert("RGB"))
        ax = axes_in(j * col_w + 0.04, y0, col_w - 0.08, row1_h)
        blank(ax)
        if key != "r110_4_axis":
            ax.imshow(img, interpolation="lanczos")
            ax.set_aspect("equal", adjustable="box")
        else:
            m = meta["panels"][key]
            wx, wy = m["win"]
            ax.imshow(img, extent=(-wx / 2, wx / 2, -wy / 2, wy / 2), interpolation="lanczos")
            ax.set_xlim(-wx / 2, m["y_hi"] + 8.0)
            ax.set_ylim(-wy / 2 - 3.0, wy / 2 + 3.0)
            ax.set_aspect("equal", adjustable="box")
            for z, t in zip(m["z_planes"], ("lower", "upper")):
                ax.text(m["y_hi"] + 0.9, z, t, ha="left", va="center", fontsize=8.0, color=FS.MUTED)
            ax.annotate("", xy=(m["y_hi"], -wy / 2 - 0.9), xytext=(m["y_lo"], -wy / 2 - 0.9),
                        arrowprops=dict(arrowstyle="<->", color=FS.INK, lw=0.7, shrinkA=0, shrinkB=0))
            ax.text(0.5 * (m["y_lo"] + m["y_hi"]), -wy / 2 - 1.3, f"$W$ = {m['W_A'] / 10:.2f} nm", ha="center",
                    va="top", fontsize=8.0, color=FS.INK)
        title(j * col_w + 0.03, y0 - 0.02, s, lab)
    # ---- rows 2 and 3: top views
    rows = (("r100_3", "r100_3h", "r100_6"), ("r110_4", "r110_4h", "r110_8"))
    letters = ("def", "ghi")
    for i, keys in enumerate(rows):
        yy = row1_h + title_h + gap_h + title_h + i * (top_h + title_h + gap_h)
        for j, key in enumerate(keys):
            m = meta["panels"][key]
            img = np.asarray(Image.open(RDIR / f"{key}.png").convert("RGB"))
            ax = axes_in(j * col_w + 0.03, yy, top_w, top_h)
            blank(ax)
            ax.imshow(img, extent=(-WIN_X / 2, WIN_X / 2, -WIN_Y / 2, WIN_Y / 2), interpolation="lanczos")
            ax.set_xlim(-WIN_X / 2, WIN_X / 2 + strip)
            ax.set_ylim(-WIN_Y / 2, WIN_Y / 2)
            ax.set_aspect("equal")
            xa = WIN_X / 2 + 1.1
            ax.annotate("", xy=(xa, m["y_hi"]), xytext=(xa, m["y_lo"]),
                        arrowprops=dict(arrowstyle="<->", color=FS.INK, lw=0.7, shrinkA=0, shrinkB=0))
            for yv in (m["y_lo"], m["y_hi"]):
                ax.plot([WIN_X / 2 + 0.2, xa + 0.6], [yv, yv], color=FS.INK, lw=0.4)
            ax.text(xa + 0.7, 0, f"$W$ = {m['W_A'] / 10:.2f} nm", rotation=90, ha="left", va="center", fontsize=8.0,
                    color=FS.INK)
            title(j * col_w + 0.03, yy - 0.02, letters[i][j], m["label"])
            if j == 0:                                    # transport direction, once per row
                ax.annotate("", xy=(-WIN_X / 2 + 6.0, -WIN_Y / 2 + 0.9), xytext=(-WIN_X / 2 + 0.6, -WIN_Y / 2 + 0.9),
                            arrowprops=dict(arrowstyle="-|>", color=FS.MUTED, lw=0.7, mutation_scale=6))
                ax.text(-WIN_X / 2 + 6.5, -WIN_Y / 2 + 0.9, f"transport [{m['direction']}]", va="center",
                        fontsize=8.0, color=FS.MUTED)
    # ---- species key
    h = [Line2D([], [], ls="", marker="o", ms=ms, mfc=COL[e], mec="#666666", mew=0.4) for e, ms in
         (("Sc", 5.2), ("N", 3.7), ("H", 2.8))]
    h.append(Line2D([], [], ls="", marker="o", ms=2.8, mfc="#f2f2f2", mec="#c8c8c8", mew=0.4))
    fig.legend(h, ["Sc", "N", "H", "H of the lower plane (top views)"], loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, 0.0), fontsize=8.0, handletextpad=0.3, columnspacing=1.6)
    FS.save(fig, OUT)


if __name__ == "__main__":
    {"render": render, "assemble": assemble}[sys.argv[1]]()
