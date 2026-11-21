"""Print style shared by the manuscript figures (Springer Nature sn-jnl, options pdflatex, sn-mathphys-num).

Each figure is drawn at the width it is included at, so the font sizes set here are the printed sizes.
Text widths were measured by a test build on 2026-11-19:
  * main text (sn-jnl): \\textwidth = 372.0 pt = 5.147 in;
  * Supplementary Material (article, 11 pt, A4, 2.5 cm margins): \\textwidth = 455.24 pt = 6.299 in.
Use `width(frac)` / `width(frac, sm=True)` with the \\linewidth fraction of the \\includegraphics call. Figures are
saved without bbox_inches="tight", which would change the width; layout is done by constrained_layout.
`check(fig)` prints the figure size and the smallest font, and refuses fonts below 8 pt (Springer: 8-12 pt).
`save` writes a vector PDF (fonts embedded, Type 42) beside the 300 dpi PNG.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.text import Text  # noqa: E402

TEXTWIDTH_PT = 372.0
SM_TEXTWIDTH_PT = 455.24
PT_PER_IN = 72.27                     # TeX points
DPI = 300
MIN_PT = 8.0                          # Springer artwork guideline: lettering 8-12 pt (raised from 6.5, 2026-11-19)

# Okabe–Ito, fixed per family (as in the T1 scripts)
COL = {"100": "#0072B2", "110": "#D55E00", "110h": "#009E73", "100h": "#CC79A7"}
SKC = "#56B4E9"
INK = "#222222"
MUTED = "#666666"
GRID = "#e6e6e6"
IRDS = 1.979


def width(frac=1.0, sm=False):
    """Figure width (in) for an \\includegraphics[width=frac\\linewidth] in the main text (or the SM)."""
    return frac * (SM_TEXTWIDTH_PT if sm else TEXTWIDTH_PT) / PT_PER_IN


def use(base=7.5):
    """rcParams for print: base font `base` pt (at least MIN_PT); ticks and legends 0.5 pt smaller (never below MIN_PT)."""
    base = max(base, MIN_PT)
    small = max(base - 0.5, MIN_PT)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": base, "axes.labelsize": base, "axes.titlesize": base,
        "xtick.labelsize": small, "ytick.labelsize": small, "legend.fontsize": small,
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 0.6, "axes.edgecolor": "#333333", "axes.labelcolor": INK, "axes.labelpad": 2.5,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.45, "ytick.minor.width": 0.45,
        "xtick.major.size": 2.8, "ytick.major.size": 2.8, "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
        "xtick.major.pad": 2.0, "ytick.major.pad": 2.0, "xtick.color": "#333333", "ytick.color": "#333333",
        "lines.linewidth": 1.0, "lines.markersize": 3.0, "lines.markeredgewidth": 0.8, "patch.linewidth": 0.5,
        "legend.frameon": False, "legend.handlelength": 2.0, "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.2, "legend.borderaxespad": 0.2, "legend.labelspacing": 0.3,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.dpi": 100, "savefig.dpi": DPI, "savefig.facecolor": "white", "figure.facecolor": "white",
        "figure.constrained_layout.use": True,
        "figure.constrained_layout.h_pad": 0.02, "figure.constrained_layout.w_pad": 0.02,
        "figure.constrained_layout.hspace": 0.02, "figure.constrained_layout.wspace": 0.02,
    })


def style(ax, grid=True):
    if grid:
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def letter(ax, s, label=None, header=None, header_color=INK, pad=3):
    """Bold panel letter "(s)" as the left title, optionally followed by a short label; `header` is a centred
    column header (e.g. the ribbon family, in its colour)."""
    ax.set_title(rf"$\mathbf{{({s})}}$" + (f" {label}" if label else ""), loc="left", color=INK, pad=pad)
    if header:
        ax.set_title(header, loc="center", color=header_color, pad=pad)


def check(fig, name):
    """Print the size and the smallest visible font; refuse anything below MIN_PT."""
    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    w, h = fig.get_size_inches()
    print(f"{name}: {w:.3f} × {h:.3f} in ({w * PT_PER_IN:.1f} pt wide), fonts {min(sizes):.2f}–{max(sizes):.2f} pt")
    assert min(sizes) >= MIN_PT - 1e-9, f"{name}: font below {MIN_PT} pt"


def save(fig, path, name=None):
    check(fig, name or path.name)
    fig.savefig(path, dpi=DPI)
    pdf = path.with_suffix(".pdf")
    fig.savefig(pdf, dpi=600)            # vector; any raster content (Fig. 1 renders) embedded at up to 600 dpi
    plt.close(fig)
    print("wrote", path, "and", pdf.name)
