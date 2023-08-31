"""

SEISMIC-RNA Logo Module

========================================================================

Draw the logo for SEISMIC-RNA as a PNG file.

"""

import unittest as ut

from functools import cache
from logging import getLogger
from pathlib import Path

from matplotlib import colormaps, pyplot as plt
import numpy as np

logger = getLogger(__name__)

LOGO_PAR_DIR_NAME = "seismic-rna"
LOGO_PAR_DIR = Path(__file__).parent
LOGO_FILE = LOGO_PAR_DIR.joinpath("logo.png")

SEGMENTS_PER_SIDE = 5040  # = 7!
SEGMENTS_FOR_INTERPOLATION = 362880  # = 9!
MAX_X = 3.

MIN_WIDTH = 3.
MAX_WIDTH = 18.

COLOR_MAP = "inferno"

TAU = 2. * np.pi


def calc_y(x: np.ndarray):
    """
    Sine times the PDF of the standard normal distribution:

    y(x) = 2π sin(2π x) N(x; 0, 1)
         = 2π sin(2π x) [exp(-x² / 2) / √(2π)]
         = √(2π) exp(-x² / 2) sin(2π x)
    """
    return np.sqrt(TAU) * np.exp(-np.square(x) / 2.) * np.sin(TAU * x)


def calc_ds(x: np.ndarray):
    """
    Length of each segment connecting each pair of consecutive x.
    """
    return np.sqrt(np.square(np.diff(x)) + np.square(np.diff(calc_y(x))))


def calc_cum_dist(x: np.ndarray):
    """
    Cumulative distance traveled on the curve at each x coordinate.
    """
    return np.concatenate([[0.], np.cumsum(calc_ds(x))])


@cache
def calc_x_interval(x_min: float, x_max: float, n_segs: int, n_interp: int):
    """
    Calculate the x values on the interval [x_min, x_max].
    """
    # Compute uniformly spaced x values for interpolating the distance.
    x_uniform = np.linspace(x_min, x_max, n_interp + 1)
    # Compute the cumulative distance at each uniform x coordinate.
    cum_dist = calc_cum_dist(x_uniform)
    # Compute uniformly spaced cumulative distances, each pair of which
    # represents one segment that will be output.
    cum_dist_uniform = np.linspace(0., cum_dist[-1], n_segs + 1)
    # Interpolate the x coordinates that correspond to these uniformly
    # spaced cumulative distances.
    return np.interp(cum_dist_uniform, cum_dist, x_uniform)


@cache
def calc_x(x_max: float, n_segs: int, n_interp: int):
    """
    Calculate the x values on the interval [-x_max, x_max].
    """
    x_half = calc_x_interval(0., x_max, n_segs, n_interp)
    return np.concatenate([-x_half[-1: -x_half.size: -1], x_half])


@cache
def points():
    """ Points of the curve. """
    x = calc_x(MAX_X, SEGMENTS_PER_SIDE, SEGMENTS_FOR_INTERPOLATION)
    return x, calc_y(x)


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
    return list(cmap(np.linspace(0., 1., len(segments()))))


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


class TestConstants(ut.TestCase):
    """ Test constants of the logo module. """

    def test_positive_segments(self):
        """ Test that SEGMENTS_PER_SIDE is a positive integer. """
        self.assertIsInstance(SEGMENTS_PER_SIDE, int)
        self.assertGreater(SEGMENTS_PER_SIDE, 0)


class TestPoints(ut.TestCase):
    """ Test function `logo.points`. """

    def test_length(self):
        """ Test that the number of points is correct. """
        xt, yt = points()
        self.assertEqual(len(xt), 2 * SEGMENTS_PER_SIDE + 1)
        self.assertEqual(len(yt), 2 * SEGMENTS_PER_SIDE + 1)

    def test_symmetric_x(self):
        """ Test that the x values are symmetric about x = 0. """
        xt = points()[0]
        self.assertTrue(np.allclose(xt, -xt[::-1]))

    def test_symmetric_y(self):
        """ Test that the y values are symmetric about x = 0. """
        yt = points()[1]
        self.assertTrue(np.allclose(yt, -yt[::-1]))


class TestSegments(ut.TestCase):
    """ Test function `logo.segments`. """

    def test_length(self):
        """ Test that segments has the correct length. """
        self.assertEqual(len(segments()), 2 * SEGMENTS_PER_SIDE)

    def test_compare_to_points(self):
        """ Test that segments accurately lists the points. """
        xt, yt = points()
        for i, ((x0, x1), (y0, y1)) in enumerate(segments()):
            self.assertEqual(x0, xt[i])
            self.assertEqual(x1, xt[i + 1])
            self.assertEqual(y0, yt[i])
            self.assertEqual(y1, yt[i + 1])

    def test_equal_lengths(self):
        """ Test that every segment has the same length. """
        # Compute each segment's length using the Pythagorean theorem.
        xt, yt = points()
        lengths = np.sqrt((xt[1:] - xt[:-1]) ** 2 + (yt[1:] - yt[:-1]) ** 2)
        # Ensure that all segments have nearly the same length.
        self.assertTrue(np.allclose(lengths, np.median(lengths),
                                    atol=1.e-6, rtol=1.e-2))


class TestWidths(ut.TestCase):
    """ Test function `logo.widths`. """

    def test_length(self):
        """ Test that the number of widths is correct. """
        self.assertEqual(widths().shape, (2 * SEGMENTS_PER_SIDE,))


class TestColors(ut.TestCase):
    """ Test function `logo.colors`. """

    def test_length(self):
        """ Test that the number of colors is correct. """
        self.assertEqual(len(colors()), 2 * SEGMENTS_PER_SIDE)


class TestDraw(ut.TestCase):
    """ Test function `logo.draw`. """

    def test_draw(self):
        """ Run draw(). """
        self.assertIsNone(draw())


if __name__ == "__main__":
    ut.main(verbosity=2)
