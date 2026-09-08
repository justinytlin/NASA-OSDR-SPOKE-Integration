"""Shared figure style (dataviz reference palette; light surface)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#ffffff"  # plain white background (matplotlib default)
GRAY_MARK = "#b3b1aa"
# fixed categorical order: one slot per tissue, never re-assigned
TISSUE_COLOR = {"hippocampus": "#2a78d6", "thymus": "#eb6834", "spleen": "#1baf7a", "liver": "#eda100",
                "kidney": "#008300", "heart": "#e34948", "adrenal gland": "#4a3aa7", "eye": "#e87ba4",
                "retina": "#e87ba4", "brain": "#2a78d6"}
TISSUE_MARKER = {"hippocampus": "o", "thymus": "s", "spleen": "^", "liver": "D", "kidney": "o", "heart": "s",
                 "adrenal gland": "^", "eye": "D", "retina": "o", "brain": "s"}
SEQ_BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
            "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
DIV_MID = "#f0efec"
DIV_BLUE, DIV_RED = "#2a78d6", "#e34948"
DATASETS = [  # display order
    ("OSD-564", "hippocampus"), ("OSD-244", "thymus"), ("OSD-289", "thymus"),
    ("OSD-246", "spleen"), ("OSD-288", "spleen"), ("OSD-245", "liver")]


def setup():
    plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
                         "font.size": 10, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                         "axes.labelcolor": INK2, "savefig.facecolor": SURFACE, "figure.facecolor": SURFACE,
                         "axes.facecolor": SURFACE})


def seq_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("seqblue", [SURFACE] + SEQ_BLUE)


def div_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("divbr", [DIV_BLUE, DIV_MID, DIV_RED])


def despine(ax, keep=("bottom",)):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)


def ptxt(p):
    return "p < 1e-15" if p < 1e-15 else f"p = {p:.2g}"
