"""

SEISMIC-RNA Logo Core Module

========================================================================

Draw the logo for SEISMIC-RNA as a PNG file.

------------------------------------------------------------------------

©2023, the Rouskin Lab.

This file is part of SEISMIC-RNA.

SEISMIC-RNA is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

SEISMIC-RNA is distributed in the hope that it will be useful, but WITH
NO WARRANTY; not even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along
with SEISMIC-RNA. If not, see https://www.gnu.org/licenses/.

========================================================================

"""

from functools import cache
from logging import getLogger
from pathlib import Path

from matplotlib import colormaps, pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

logger = getLogger(__name__)

LOGO_PAR_DIR_NAME = "seismic-rna"
LOGO_PAR_DIR = Path(__file__).parent.parent.parent.parent
LOGO_FILE = LOGO_PAR_DIR.joinpath("logo.png")

SEGMENTS = 1080

MIN_WIDTH = 3.
MAX_WIDTH = 18.

COLOR_MAP = "inferno"


@cache
def points():
    """ Return the points of the curve. """

    tau = 2. * np.pi

    def y_of_x(x: np.ndarray):
        """
        Sine times the PDF of the standard normal distribution:

        y(x) = 2π sin(2π x) N(x; 0, 1)
             = 2π sin(2π x) [exp(-x² / 2) / √(2π)]
             = √(2π) sin(2π x) exp(-x² / 2)
        """
        return np.sqrt(tau) * np.sin(tau * x) * np.exp(-x ** 2 / 2.)

    def x_of_t(t_eval: np.ndarray):
        """ Integrate dx/dt numerically to find x(t). """

        def dx_dt(_: float, x: np.ndarray):
            """
            Compute the distance traveled along the curve per distance
            traveled along the x-axis.

            dy/dx = √(2π) exp(-x² / 2) [2π cos(2π x) - x sin(2π x)]
            ds/dx = √[(dy/dx)² + (dx/dx)²]
                  = √(2π exp(-x²) [2π cos(2π x) - x sin(2π x)]² + 1)

            Set the speed (time derivative of distance along the curve)
            to a constant and solve for the time derivative of x.

            ds/dt = c  (constant)
            dx/dt = (ds/dt) / (ds/dx)
                  = c / √(2π exp(-x²) [2π cos(2π x) - x sin(2π x)]² + 1)
            """
            # Arbitrary constant: value is irrelevant so long as c > 0.
            c = 1.
            # Evaluate the inner part of the function.
            tau_x = tau * x
            a = tau * np.cos(tau_x) - x * np.sin(tau_x)
            # Compute the time derivative of x as a function of x.
            return c / np.sqrt((tau * np.exp(-x ** 2) * a ** 2) + 1.)

        # Integrate the differential equation numerically.
        t_span = t_eval[0], t_eval[-1]
        solution = solve_ivp(dx_dt, t_span, [0.], t_eval=t_eval)
        x_eval = solution.y[0]
        return x_eval

    t = np.linspace(0., 13., SEGMENTS // 2 + 1)
    xt = x_of_t(t)
    xt = np.hstack([-xt[-1: -xt.size: -1], xt])
    yt = y_of_x(xt)
    return xt, yt


def _segments():
    """ Yield each line segment. """
    xs, ys = points()
    x0, y0 = xs[0], ys[0]
    for x1, y1 in zip(xs[1:], ys[1:], strict=True):
        yield (float(x0), float(x1)), (float(y0), float(y1))
        x0, y0 = x1, y1


@cache
def segments():
    """ Return a list of each line segment. """
    return list(_segments())


@cache
def widths():
    """ Return the width of each segment. """
    # Compute the linear component of the line width.
    xs = points()[0]
    mean_xs = (xs[1:] + xs[:-1]) / 2.
    linear = np.exp(-mean_xs ** 2 / 2.)
    # Normalize the linear component to the range [0, 1].
    norm_linear = (linear - np.min(linear)) / (np.max(linear) - np.min(linear))
    # Transform the linear component based on the desired dimensions.
    return (MAX_WIDTH - MIN_WIDTH) * norm_linear + MIN_WIDTH


@cache
def colors():
    """ Return the color of each segment. """
    cmap = colormaps[COLOR_MAP]
    return list(cmap(np.linspace(0., 1., SEGMENTS)))


def draw():
    """ Draw the logo for SEISMIC-RNA. """
    if LOGO_PAR_DIR.name != LOGO_PAR_DIR_NAME:
        logger.error(f"The logo must be written into the directory named "
                     f"{LOGO_PAR_DIR_NAME}, but the directory in this file "
                     f"system is named {LOGO_PAR_DIR.name}. Most likely, you "
                     f"are trying to generate the logo from an installed copy "
                     f"of SEISMIC-RNA. The logo can only be generated from the "
                     f"source code before it is installed.")
        return
    fig, ax = plt.subplots()
    for seg, w, c in zip(segments(), widths(), colors(),
                         strict=True):
        ax.plot(*seg, color=c, linewidth=w, solid_capstyle="round")
    ax.set_aspect(1.0)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(LOGO_FILE, dpi=300, transparent=True)
    plt.close()


if __name__ == "__main__":
    draw()
