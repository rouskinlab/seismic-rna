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
from scipy.integrate import solve_ivp

LOGO_PAR_DIR = Path(__file__).parent.parent.parent.parent
LOGO_FILE = LOGO_PAR_DIR.joinpath("logo.png")

SEGMENTS = 1080

MIN_WIDTH = 3.
MAX_WIDTH = 18.


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
        return np.sqrt(tau) * np.sin(tau * x) * np.exp(-x**2 / 2)

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
            return c / np.sqrt((tau * np.exp(-x * x) * a * a) + 1.)

        # Integrate the differential equation numerically.
        t_span = t_eval[0], t_eval[-1]
        solution = solve_ivp(dx_dt, t_span, [0.], t_eval=t_eval)
        x_eval = solution.y[0]
        return x_eval

    t = np.linspace(0., 4. * np.pi, SEGMENTS // 2 + 1)
    xt = x_of_t(t)
    xt = np.hstack([-xt[-1: -xt.size: -1], xt])
    yt = y_of_x(xt)
    return xt, yt


def segments():
    """ Return an iterator for each line segment. """
    xs, ys = points()
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
    return list(cmap(np.linspace(0., 1., SEGMENTS)))


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
