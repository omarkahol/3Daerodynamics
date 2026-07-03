import matplotlib.pyplot as plt

def setup_plot_style():
    """Sets professional plotting style configurations for Matplotlib."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.figsize": (5.5, 4.0),
        "lines.linewidth": 1.5,
        "savefig.bbox": "tight",
        "savefig.format": "pdf"
    })
