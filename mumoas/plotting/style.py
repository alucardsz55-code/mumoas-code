from __future__ import annotations


def set_publication_style(font_scale: float = 1.0) -> None:
    """Apply restrained defaults suitable for manuscript-review figures."""
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 9 * font_scale,
            "axes.titlesize": 10 * font_scale,
            "figure.dpi": 150,
            "font.size": 9 * font_scale,
            "legend.frameon": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 300,
            "xtick.labelsize": 8 * font_scale,
            "ytick.labelsize": 8 * font_scale,
        }
    )
