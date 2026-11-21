"""OVITO renders of bulk rocksalt ScN and the M-ScN sheet (precursor of plan §6 Fig. 1: structures and W).

Same recipe as the Cr2N and ns projects: Tachyon ray tracer with ambient occlusion, orthographic camera, panels
assembled with matplotlib.

    LD_LIBRARY_PATH=~/.venvs/structviz/syslib/root/usr/lib/x86_64-linux-gnu \\
        ~/.venvs/structviz/bin/python scripts/render_structures.py [--draft]

Geometries (ideal, built here; no DFT output is read yet):
  * bulk rocksalt ScN, a = 4.50 Å (experiment 4.501 Å; the HSE06 references use 4.48–4.50);
  * M-ScN, the two-plane rocksalt-(001) bilayer (P4/nmm), primitive a = 3.11 Å (Sc–N in plane 2.20 Å, Liang's
    4-f.u. cell 4.40 Å) and interlayer distance 2.18 Å: Liang 2021 PBE-D3, the same values as the A0 sheet inputs.
    Our relaxed A2 geometry replaces them when it exists;
  * the 3-plane slab (hold-out A3), ideal stacking with the same spacings: illustrative, unrelaxed.
Writes figures/Fig1_structures_bulk_sheet.png.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.transforms import ScaledTranslation
from PIL import Image
from ase import Atoms

from ovito.io.ase import ase_to_ovito
from ovito.modifiers import CreateBondsModifier
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import BondsVis, Viewport, TachyonRenderer

DRAFT = "--draft" in sys.argv
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "Fig1_structures_bulk_sheet.png"
TMP = Path(tempfile.mkdtemp(prefix="ovito_"))

A_BULK = 4.50                  # Å
A_PRIM = 3.11                  # Å, M-ScN primitive (P4/nmm)
A_CONV = A_PRIM * np.sqrt(2)   # 4.398 Å, Liang's 4-f.u. cell
D_LAYER = 2.18                 # Å, interlayer distance (Liang PBE-D3)

COL = {"Sc": "#c9b48c", "N": "#4a64d6"}      # sand Sc, blue N: distinct in hue and lightness
RAD = {"Sc": 0.56, "N": 0.40}                # display radii, Å
INK, MUTED = "#1d1f22", "#6b6e73"
SCALE = 2 if DRAFT else 4                     # render pixels per figure pixel step


def renderer():
    if DRAFT:
        return TachyonRenderer(ambient_occlusion=False, antialiasing_samples=4, shadows=False)
    return TachyonRenderer(ambient_occlusion=True, ambient_occlusion_brightness=0.85, ambient_occlusion_samples=16,
                           antialiasing_samples=9, shadows=False, direct_light_intensity=0.9)


def pipeline(pos, sym, cell=None):
    """Static pipeline. cell: optional (3, 3) rows drawn as an outline at the origin."""
    big = np.diag([80.0, 80.0, 80.0])
    dc = ase_to_ovito(Atoms(symbols=sym, positions=pos, cell=big if cell is None else cell, pbc=False))
    if cell is None:
        dc.cell_.vis.enabled = False
    else:
        dc.cell_.vis.line_width = 0.07
        dc.cell_.vis.rendering_color = to_rgb("#55585e")
    types = dc.particles_.particle_types_
    for el in set(sym):
        t = types.type_by_name_(el)
        t.color, t.radius = to_rgb(COL[el]), RAD[el]
    pipe = Pipeline(source=StaticSource(data=dc))
    b = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
    b.set_pairwise_cutoff("Sc", "N", 2.6)          # nearest neighbours only (2.18–2.25 Å)
    b.vis.width = 0.15
    b.vis.coloring_mode = BondsVis.ColoringMode.Uniform
    b.vis.color = to_rgb("#9a9ca1")
    pipe.modifiers.append(b)
    return pipe


def snapshot(pos, sym, direction, up, center, half_h, size, cell=None):
    """Orthographic render; scale = size[1] / (2 half_h) px per Å (keep size[0] ≥ size[1])."""
    pipe = pipeline(pos, sym, cell)
    pipe.add_to_scene()
    d = np.asarray(direction, float)
    d /= np.linalg.norm(d)
    vp = Viewport(type=Viewport.Type.Ortho, camera_dir=tuple(d),
                  camera_pos=tuple(np.asarray(center, float) - 80 * d), fov=half_h)
    vp.camera_up = tuple(up)
    f = TMP / "r.png"
    vp.render_image(size=tuple(int(v) for v in size), filename=str(f), background=(1, 1, 1), renderer=renderer())
    pipe.remove_from_scene()
    return np.asarray(Image.open(f).convert("RGB"))


# ---------------------------------------------------------------- geometries
def bulk_cube(a=A_BULK):
    """Conventional rocksalt cube with all face, edge and corner sites (27 atoms; Sc at the corners)."""
    pos, sym = [], []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                pos.append(np.array([i, j, k]) * a / 2)
                sym.append("Sc" if (i + j + k) % 2 == 0 else "N")
    return np.array(pos), sym


def sheet(nplanes=2, n=8, a=A_CONV, d=D_LAYER):
    """Rocksalt-(001) slab of `nplanes` planes on an (n+1)×(n+1) grid of spacing a/2 (a = 4-f.u. cell edge).
    Plane p swaps species relative to plane p − 1 (rocksalt stacking)."""
    pos, sym = [], []
    for p in range(nplanes):
        for i in range(n + 1):
            for j in range(n + 1):
                pos.append([i * a / 2, j * a / 2, p * d])
                sym.append("Sc" if (i + j + p) % 2 == 0 else "N")
    return np.array(pos), sym


# ---------------------------------------------------------------- renders
def main():
    panels = {}
    # (a) bulk cube, oblique
    pos, sym = bulk_cube()
    c = np.full(3, A_BULK / 2)
    panels["a"] = snapshot(pos, sym, (-1.0, -1.45, -0.75), (0, 0, 1), c, 4.35, (380 * SCALE, 360 * SCALE),
                           cell=np.eye(3) * A_BULK)
    # (b) M-ScN sheet, oblique (a 4×4 patch of the 4.40 Å cell)
    pos, sym = sheet(2, 5)
    c = [pos[:, 0].mean(), pos[:, 1].mean(), D_LAYER / 2]
    panels["b"] = snapshot(pos, sym, (-1.0, -1.55, -0.55), (0, 0, 1), c, 8.1, (380 * SCALE, 360 * SCALE))
    # (c) top view of the sheet
    pos_t, sym_t = sheet(2, 6)
    top = pos_t[:, 2] > D_LAYER / 2                        # the lower plane hides exactly under it, species swapped
    pos_t, sym_t = pos_t[top], [x for x, k in zip(sym_t, top) if k]
    ct = [pos_t[:, 0].mean(), pos_t[:, 1].mean(), D_LAYER]
    half_c = 8.3
    panels["c"] = snapshot(pos_t, sym_t, (0, 0, -1), (0, 1, 0), ct, half_c, (380 * SCALE, 360 * SCALE))
    # (d, e) side views along the conventional [110] = primitive a1: each column is a single species
    side = []
    for npl in (2, 3):
        pos_s, sym_s = sheet(npl, 12)
        keep = pos_s[:, 0] - pos_s[:, 1] >= -1e-6           # half the patch: columns from one edge
        cs = [0, 0, (npl - 1) * D_LAYER / 2]
        u = np.array([1.0, -1.0, 0]) / np.sqrt(2)
        proj = pos_s[keep] @ u
        cs = u * (proj.mean() + 4.4) + np.array([0, 0, (npl - 1) * D_LAYER / 2])
        side.append((snapshot(pos_s[keep], [s for s, k in zip(sym_s, keep) if k], (1, 1, 0), (0, 0, 1), cs, 3.9,
                              (760 * SCALE, 210 * SCALE)), proj, npl))

    # ------------------------------------------------------------ assemble
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.5], hspace=0.2, wspace=0.08,
                          left=0.02, right=0.98, top=0.95, bottom=0.03)

    def show(ax, img, label, title):
        ax.imshow(img, interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.text(0.0, 1.02, label, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="bottom")
        ax.text(0.0, 1.02, title, transform=ax.transAxes + ScaledTranslation(25 / 72, 0, fig.dpi_scale_trans),
                fontsize=8.5, color=INK, va="bottom")

    ax = fig.add_subplot(gs[0, 0])
    show(ax, panels["a"], "(a)", f"bulk rocksalt ScN, Fm3̄m, a = {A_BULK:.2f} Å")
    ax = fig.add_subplot(gs[0, 1])
    show(ax, panels["b"], "(b)", "M-ScN sheet: rocksalt-(001) bilayer, P4/nmm")

    # (c) with cell outlines
    ax = fig.add_subplot(gs[1, 0])
    img = panels["c"]
    show(ax, img, "(c)", "M-ScN, top view (upper plane)")
    H, W = img.shape[:2]
    s = H / (2 * half_c)

    def px(x, y):
        return W / 2 + (x - ct[0]) * s, H / 2 - (y - ct[1]) * s
    o = np.array([2 * A_CONV / 2 + A_CONV / 2 * 0, 2 * A_CONV / 2])       # a Sc site near the centre-left
    o = np.array([ct[0] - A_CONV / 2, ct[1] - A_CONV / 2])
    a1, a2 = np.array([A_CONV / 2, A_CONV / 2]), np.array([-A_CONV / 2, A_CONV / 2])     # primitive, 3.11 Å
    prim = np.array([o, o + a1, o + a1 + a2, o + a2, o]) + np.array([A_CONV / 2, 0])
    conv = np.array([o, o + [A_CONV, 0], o + [A_CONV, A_CONV], o + [0, A_CONV], o])
    P = np.array([px(*p) for p in prim]); Cv = np.array([px(*p) for p in conv])
    ax.plot(Cv[:, 0], Cv[:, 1], color=INK, lw=1.1, ls=(0, (4, 2.5)))
    ax.plot(P[:, 0], P[:, 1], color=INK, lw=1.6)
    ax.set_xlim(0, W); ax.set_ylim(H, 0)

    # legend panel next to (c)
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    yy = 0.86
    for el, lab in (("Sc", "Sc"), ("N", "N")):
        ax.scatter([0.08], [yy], s=260 if el == "Sc" else 150, color=COL[el], edgecolor="#55585e", lw=0.6)
        ax.text(0.16, yy, lab, va="center", fontsize=9, color=INK)
        yy -= 0.12
    notes = [("solid square", f"primitive cell, a = {A_PRIM:.2f} Å (2 f.u.)"),
             ("dashed square", f"Liang's 4-f.u. cell, {A_CONV:.2f} Å"),
             ("in-plane Sc–N", f"{A_PRIM / np.sqrt(2):.2f} Å (bulk {A_BULK / 2:.3f} Å)"),
             ("interlayer Sc–N", f"{D_LAYER:.2f} Å"),
             ("lower plane", "directly beneath (c), Sc and N swapped"),
             ("geometry", "Liang 2021 PBE-D3 (= A0 inputs);\nour A2 relaxation replaces it")]
    yy -= 0.02
    for k, v in notes:
        ax.text(0.02, yy, k, va="top", fontsize=7.5, color=MUTED)
        ax.text(0.36, yy, v, va="top", fontsize=7.5, color=INK)
        yy -= 0.12 if "\n" not in v else 0.16
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # (d), (e) side views with the interlayer spacing marked
    sub = gs[2, :].subgridspec(1, 2, wspace=0.06)
    for k, (img, proj, npl) in enumerate(side):
        ax = fig.add_subplot(sub[0, k])
        show(ax, img, "(d)" if k == 0 else "(e)",
             "M-ScN, side view along [110]" if npl == 2 else "3-plane slab (hold-out A3), ideal stacking")
        H, W = img.shape[:2]
        s = H / (2 * 3.9)
        zc = (npl - 1) * D_LAYER / 2
        xr = W - 0.06 * W
        for p in range(npl - 1):
            y0, y1 = H / 2 - (p * D_LAYER - zc) * s, H / 2 - ((p + 1) * D_LAYER - zc) * s
            ax.annotate("", xy=(xr, y1), xytext=(xr, y0),
                        arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9, shrinkA=0, shrinkB=0))
            ax.text(xr - 0.012 * W, (y0 + y1) / 2, f"{D_LAYER:.2f} Å", ha="right", va="center", fontsize=7.5,
                    color=INK, bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
        ax.set_xlim(0, W); ax.set_ylim(H, 0)

    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=200 if DRAFT else 300, facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
