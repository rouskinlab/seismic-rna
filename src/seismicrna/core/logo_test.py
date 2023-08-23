"""

Tests for SEISMIC-RNA Logo Core Module

========================================================================

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

import unittest as ut

import numpy as np

from .logo import SEGMENTS, points, segments, widths, colors


class TestConstants(ut.TestCase):
    """ Test constants of the logo module. """

    def test_positive_even_segments(self):
        """ Test that SEGMENTS is an even positive integer. """
        self.assertIsInstance(SEGMENTS, int)
        self.assertGreater(SEGMENTS, 0)
        self.assertEqual(SEGMENTS % 2, 0)


class TestPoints(ut.TestCase):
    """ Test function `logo.points`. """

    def test_length(self):
        """ Test that the number of points is correct. """
        xt, yt = points()
        self.assertEqual(len(xt), SEGMENTS + 1)
        self.assertEqual(len(yt), SEGMENTS + 1)

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
        self.assertEqual(len(segments()), SEGMENTS)

    def test_compare_to_points(self):
        """ Test that segments accurately lists the points. """
        xt, yt = points()
        for i, ((x0, x1), (y0, y1)) in enumerate(segments()):
            self.assertEqual(x0, xt[i])
            self.assertEqual(x1, xt[i + 1])
            self.assertEqual(y0, yt[i])
            self.assertEqual(y1, yt[i + 1])

    '''
    def test_equal_lengths(self):
        """ Test that every segment has the same Euclidian length. """
        # Compute each segment's length using the Pythagorean theorem.
        xt, yt = points()
        lengths = np.sqrt((xt[1:] - xt[:-1]) ** 2 + (yt[1:] - yt[:-1]) ** 2)
        # Ensure that all lengths are practically equal.
        self.assertTrue(np.allclose(lengths, np.median(lengths)))
    '''


class TestWidths(ut.TestCase):
    """ Test function `logo.widths`. """

    def test_length(self):
        """ Test that the number of widths is correct. """
        self.assertEqual(widths().shape, (SEGMENTS,))


class TestColors(ut.TestCase):
    """ Test function `logo.widths`. """

    def test_length(self):
        """ Test that the number of colors is correct. """
        self.assertEqual(len(colors()), SEGMENTS)
