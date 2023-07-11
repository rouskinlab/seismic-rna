"""
Logo Module
========================================================================
Auth: Matty

Draw the logo for SEISMIC-RNA.
"""

from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np

LOGO_PAR_DIR = Path(__file__).parent.parent.parent.parent
LOGO_FILE = LOGO_PAR_DIR.joinpath("logo.png")

SEGMENTS = 1080

MIN_WIDTH = 3.
MAX_WIDTH = 18.


def ft():
    """
    Return sine times the PDF of the standard normal distribution:

    f(x) = 2π sin(2πx) N(x; 0, 1) = 2π sin(2πx) exp(-x^2 / 2) / √2π
         = √2π sin(2πx) exp(-x^2 / 2)

    Convert to a parametric distribution so that each segment can be the
    same length along the curve, not the same width along the x-axis.

    

    """
    tau = 2. * np.pi
    x = np.linspace(-3., 3., SEGMENTS + 1)
    y = np.sqrt(tau) * np.sin(tau * x) * np.exp(-0.5 * x * x)
    return x, y


def segments():
    """ Return an iterator for each line segment. """
    xs, ys = fx()
    x0, y0 = xs[0], ys[0]
    for x1, y1 in zip(xs[1:], ys[1:]):
        yield [x0, x1], [y0, y1]
        x0, y0 = x1, y1


def widths():
    """ Return the width of each segment. """
    return np.hstack([np.linspace(MIN_WIDTH, MAX_WIDTH, SEGMENTS // 2),
                      np.linspace(MAX_WIDTH, MIN_WIDTH, SEGMENTS // 2)])


def colors():
    """ Return the color of each segment. """
    cmap = get_cmap("inferno")
    return list(cmap(np.hstack([np.linspace(1., 0., SEGMENTS // 2),
                                np.linspace(0., 1., SEGMENTS // 2)])))


def draw():
    """ Draw the logo for SEISMIC-RNA. """
    fig, ax = plt.subplots()
    for seg, w, c in zip(segments(), widths(), colors(),
                         strict=True):
        ax.plot(*seg, color=c, linewidth=w, solid_capstyle="round")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(LOGO_FILE, dpi=300, transparent=True)
    plt.close()


if __name__ == "__main__":
    draw()
