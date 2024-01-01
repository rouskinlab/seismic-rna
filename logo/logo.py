"""

SEISMIC-RNA Logo Module

========================================================================

Draw the logo for SEISMIC-RNA as a PNG file.

"""

import unittest as ut

from functools import cache
from logging import getLogger

import numpy as np
from matplotlib import colormaps, patches, pyplot as plt
from PIL import Image

logger = getLogger(__name__)

BASE_LOGO = "logo-base.png"
BASE_BLUE = "logo-blue.png"
BASE_DIMENSION = 4800
LOGO_TEMPLATE = "logo-{}.png"
LOGO_EDGES = [1200, 200]
FAVICON_TEMPLATE = "favicon-{}.ico"
FAVICON_SIZES = [(32, 32)]
BG_COLOR = "#2e8ece"

SEGMENTS_PER_SIDE = 5040  # = 7!
SEGMENTS_FOR_INTERPOLATION = 362880  # = 9!
MAX_X = 3.

MIN_WIDTH = 3.
MAX_WIDTH = 18.
MARGIN = 0.2

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


def draw_base_logo():
    """ Draw the logo for SEISMIC-RNA. """
    fig, ax = plt.subplots()
    # Plot each segment with its own line width and color.
    for seg, w, c in zip(segments(), widths(), colors(),
                         strict=True):
        ax.plot(*seg, color=c, linewidth=w, solid_capstyle="round")
    # Calculate the width and height of the graph and set the length of
    # each side of the square.
    width_inches, height_inches = (np.max(xy) - np.min(xy) for xy in points())
    length = max(width_inches, height_inches) * (1. + MARGIN)
    half_length = length / 2.
    # Hide the spines and ticks.
    plt.axis("off")
    # Set the axis dimensions and aspect ratio.
    ax.set_xlim(-half_length, half_length)
    ax.set_ylim(-half_length, half_length)
    ax.set_aspect(1.)
    # Set the figure size and margins.
    fig.subplots_adjust(left=0., right=1., top=1., bottom=0.)
    fig.set_size_inches(w=length, h=length)
    # Set the image resolution.
    dpi = BASE_DIMENSION / length
    # Save the figure.
    plt.savefig(BASE_LOGO, dpi=dpi, transparent=True)
    # Add a blue circle and save the figure again.
    ax.add_patch(patches.Circle((0., 0.), radius=half_length, fill=BG_COLOR))
    plt.savefig(BASE_BLUE, dpi=dpi, transparent=True)
    plt.close()


def draw_extra_logos():
    image = Image.open(BASE_LOGO)
    for edge in LOGO_EDGES:
        width, height = image.size
        resized = image.resize((edge, round(edge * height / width)),
                               resample=Image.LANCZOS)
        resized.save(LOGO_TEMPLATE.format(edge), format="PNG")


def draw_favicon():
    """ Draw the favicon for SEISMIC-RNA. """
    image = Image.open(BASE_BLUE)
    for size in FAVICON_SIZES:
        image.save(FAVICON_TEMPLATE.format('x'.join(map(str, size))),
                   format="ICO", sizes=[size])


def draw():
    draw_base_logo()
    draw_extra_logos()
    draw_favicon()


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

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
