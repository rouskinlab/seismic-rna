import unittest as ut
from itertools import product

import numpy as np
import pandas as pd

from seismicrna.core.mu import reframe, reframe_like

rng = np.random.default_rng()


def broadcastable(vshape: tuple[int, ...], tshape: tuple[int, ...]):
    """ Check whether values in the shape of `vshape` can be broadcast
    to the target shape of `tshape`.

    Parameters
    ----------
    vshape: tuple[int, ...]
        Shape of the values to broadcast.
    tshape: tuple[int, ...]
        Shape of the target to which the values should be broadcast.

    Returns
    -------
    bool
        Whether the value shape can be broadcast to the target shape.
    """
    if len(tshape) < len(vshape):
        # An array cannot be broadcast to a shape with fewer dimensions.
        return False
    if len(tshape) > len(vshape):
        # If the target shape has more dimensions than the array, then
        # prepend extra dimensions (each of length 1) to the array.
        return broadcastable((len(tshape) - len(vshape)) * (1,) + vshape,
                             tshape)
    # The array is broadcastable to the target shape if and only if, for
    # every dimension, the length of the original axis is 1 or equals
    # the length of the target axis.
    for vlen, tlen in zip(vshape, tshape, strict=True):
        if vlen != 1 and vlen != tlen:
            return False
    return True


class TestReframe(ut.TestCase):

    def test_float_none(self):
        for value in np.linspace(0., 1., 3):
            frame = reframe(value)
            self.assertIsInstance(frame, np.ndarray)
            self.assertEqual(frame.ndim, 0)
            self.assertTrue(np.array_equal(frame, np.array(value)))

    def test_float_ints(self):
        for ndim in range(4):
            for shape in product(range(4), repeat=ndim):
                self.assertIsInstance(shape, tuple)
                self.assertEqual(len(shape), ndim)
                for value in np.linspace(0., 1., 3):
                    frame = reframe(value, shape)
                    self.assertIsInstance(frame, np.ndarray)
                    self.assertEqual(frame.shape, shape)
                    self.assertTrue(np.allclose(frame, value))

    def test_float_index(self):
        for length in range(5):
            index = rng.integers(10, size=length)
            self.assertEqual(index.shape, (length,))
            for value in np.linspace(0., 1., 3):
                frame = reframe(value, (index,))
                self.assertIsInstance(frame, pd.Series)
                self.assertEqual(frame.shape, (length,))
                self.assertTrue(np.allclose(frame, value))
                self.assertIsInstance(frame.index, pd.Index)
                self.assertTrue(np.all(frame.index == index))

    def test_float_index_int(self):
        for nrow in range(5):
            rows = rng.integers(10, size=nrow)
            self.assertEqual(rows.shape, (nrow,))
            for ncol in range(3):
                for value in np.linspace(0., 1., 3):
                    frame = reframe(value, (rows, ncol))
                    self.assertIsInstance(frame, pd.DataFrame)
                    self.assertEqual(frame.shape, (nrow, ncol))
                    self.assertTrue(np.allclose(frame, value))
                    self.assertIsInstance(frame.index, pd.Index)
                    self.assertTrue(np.all(frame.index == rows))
                    self.assertIsInstance(frame.columns, pd.RangeIndex)
                    self.assertTrue(np.all(frame.columns == np.arange(ncol)))

    def test_float_int_index(self):
        for nrow in range(5):
            for ncol in range(3):
                cols = rng.integers(10, size=ncol)
                self.assertEqual(cols.shape, (ncol,))
                for value in np.linspace(0., 1., 3):
                    frame = reframe(value, (nrow, cols))
                    self.assertIsInstance(frame, pd.DataFrame)
                    self.assertEqual(frame.shape, (nrow, ncol))
                    self.assertTrue(np.allclose(frame, value))
                    self.assertIsInstance(frame.index, pd.RangeIndex)
                    self.assertTrue(np.all(frame.index == np.arange(nrow)))
                    self.assertIsInstance(frame.columns, pd.Index)
                    self.assertTrue(np.all(frame.columns == cols))

    def test_float_index_index(self):
        for nrow in range(5):
            rows = rng.integers(10, size=nrow)
            self.assertEqual(rows.shape, (nrow,))
            for ncol in range(3):
                cols = rng.integers(10, size=ncol)
                self.assertEqual(cols.shape, (ncol,))
                for value in np.linspace(0., 1., 3):
                    frame = reframe(value, (rows, cols))
                    self.assertIsInstance(frame, pd.DataFrame)
                    self.assertEqual(frame.shape, (nrow, ncol))
                    self.assertTrue(np.allclose(frame, value))
                    self.assertIsInstance(frame.index, pd.Index)
                    self.assertTrue(np.all(frame.index == rows))
                    self.assertIsInstance(frame.columns, pd.Index)
                    self.assertTrue(np.all(frame.columns == cols))

    def test_float_index_index_int(self):
        for nrow in range(5):
            rows = rng.integers(10, size=nrow)
            self.assertEqual(rows.shape, (nrow,))
            for ncol in range(3):
                cols = rng.integers(10, size=ncol)
                self.assertEqual(cols.shape, (ncol,))
                for nlev in range(2):
                    for value in np.linspace(0., 1., 3):
                        self.assertRaisesRegex(ValueError,
                                               "A Pandas object must have 1 "
                                               "or 2 axes, but got 3",
                                               reframe,
                                               value,
                                               (rows, cols, nlev))

    def test_float_index_index_index(self):
        for nrow in range(5):
            rows = rng.integers(10, size=nrow)
            self.assertEqual(rows.shape, (nrow,))
            for ncol in range(3):
                cols = rng.integers(10, size=ncol)
                self.assertEqual(cols.shape, (ncol,))
                for nlev in range(2):
                    levs = rng.integers(10, size=nlev)
                    self.assertEqual(levs.shape, (nlev,))
                    for value in np.linspace(0., 1., 3):
                        self.assertRaisesRegex(ValueError,
                                               "A Pandas object must have 1 "
                                               "or 2 axes, but got 3",
                                               reframe,
                                               value,
                                               (rows, cols, levs))

    def test_array_none(self):
        for ndim in range(4):
            for shape in product(range(4), repeat=ndim):
                self.assertIsInstance(shape, tuple)
                self.assertEqual(len(shape), ndim)
                value = rng.random(shape)
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, shape)
                frame = reframe(value, shape)
                self.assertIsInstance(frame, np.ndarray)
                self.assertEqual(frame.shape, shape)
                self.assertTrue(np.allclose(frame, value))

    def test_array_ints(self):
        for vdim in range(4):
            for vshape in product(range(4), repeat=vdim):
                self.assertIsInstance(vshape, tuple)
                self.assertEqual(len(vshape), vdim)
                value = rng.random(vshape)
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, vshape)
                for tdim in range(4):
                    for tshape in product(range(4), repeat=tdim):
                        self.assertIsInstance(tshape, tuple)
                        self.assertEqual(len(tshape), tdim)
                        if broadcastable(vshape, tshape):
                            fshape = np.broadcast_shapes(vshape, tshape)
                            frame = reframe(value, tshape)
                            self.assertIsInstance(frame, np.ndarray)
                            self.assertEqual(frame.shape, fshape)
                            self.assertTrue(np.allclose(frame, value))
                        else:
                            if tdim == 0:
                                err = ("cannot broadcast a non-scalar to a "
                                       "scalar array")
                            elif vdim > tdim:
                                err = ("input operand has more dimensions "
                                       "than allowed by the axis remapping")
                            else:
                                err = "could not be broadcast together"
                            self.assertRaisesRegex(ValueError,
                                                   err,
                                                   reframe,
                                                   value,
                                                   tshape)

    def test_array_index(self):
        for ndim in range(4):
            for shape in product(range(4), repeat=ndim):
                self.assertIsInstance(shape, tuple)
                self.assertEqual(len(shape), ndim)
                value = rng.random(shape)
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, shape)
                for length in range(5):
                    index = rng.integers(10, size=length)
                    self.assertEqual(index.shape, (length,))
                    if ndim == 0 or length == shape[0]:
                        if ndim <= 1:
                            frame = reframe(value, (index,))
                            self.assertIsInstance(frame, pd.Series)
                            self.assertEqual(frame.shape, (length,))
                            self.assertTrue(np.allclose(frame, value))
                            self.assertIsInstance(frame.index, pd.Index)
                            self.assertTrue(np.all(frame.index == index))
                        else:
                            self.assertRaisesRegex(ValueError,
                                                   "Data must be 1-dimensional",
                                                   reframe,
                                                   value,
                                                   (index,))
                    else:
                        self.assertRaisesRegex(ValueError,
                                               "does not match length of index",
                                               reframe,
                                               value,
                                               (index,))

    def test_array_index_int(self):
        for ndim in range(4):
            for shape in product(range(4), repeat=ndim):
                self.assertIsInstance(shape, tuple)
                self.assertEqual(len(shape), ndim)
                value = rng.random(shape)
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, shape)
                for nrow in range(5):
                    rows = rng.integers(10, size=nrow)
                    self.assertEqual(rows.shape, (nrow,))
                    for ncol in range(3):
                        if 1 <= ndim <= 2:
                            if ((ndim, nrow, ncol) == (1, shape[0], 1)
                                    or (nrow, ncol) == shape):
                                frame = reframe(value, (rows, ncol))
                                self.assertIsInstance(frame, pd.DataFrame)
                                self.assertEqual(frame.shape, (nrow, ncol))
                                if ndim == 1:
                                    self.assertTrue(np.allclose(frame.T, value))
                                else:
                                    self.assertTrue(np.allclose(frame, value))
                                self.assertIsInstance(frame.index, pd.Index)
                                self.assertTrue(np.all(frame.index == rows))
                                self.assertIsInstance(frame.columns,
                                                      pd.RangeIndex)
                                self.assertTrue(np.all(frame.columns
                                                       == np.arange(ncol)))
                            else:
                                err = ("Empty data passed with indices"
                                       if shape[0] == 0 and nrow > 0
                                       else "Shape of passed values")
                                self.assertRaisesRegex(ValueError,
                                                       err,
                                                       reframe,
                                                       value,
                                                       (rows, ncol))
                        else:
                            self.assertRaisesRegex(ValueError,
                                                   "Must pass 2-d input",
                                                   reframe,
                                                   value,
                                                   (rows, ncol))

    def test_array_int_index(self):
        for ndim in range(4):
            for shape in product(range(4), repeat=ndim):
                self.assertIsInstance(shape, tuple)
                self.assertEqual(len(shape), ndim)
                value = rng.random(shape)
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, shape)
                for nrow in range(5):
                    for ncol in range(3):
                        cols = rng.integers(10, size=ncol)
                        self.assertEqual(cols.shape, (ncol,))
                        if 1 <= ndim <= 2:
                            if ((ndim, nrow, ncol) == (1, shape[0], 1)
                                    or (nrow, ncol) == shape):
                                frame = reframe(value, (nrow, cols))
                                self.assertIsInstance(frame, pd.DataFrame)
                                self.assertEqual(frame.shape, (nrow, ncol))
                                if ndim == 1:
                                    self.assertTrue(np.allclose(frame.T, value))
                                else:
                                    self.assertTrue(np.allclose(frame, value))
                                self.assertIsInstance(frame.index,
                                                      pd.RangeIndex)
                                self.assertTrue(np.all(frame.index
                                                       == np.arange(nrow)))
                                self.assertIsInstance(frame.columns, pd.Index)
                                self.assertTrue(np.all(frame.columns == cols))
                            else:
                                err = ("Empty data passed with indices"
                                       if shape[0] == 0 and nrow > 0
                                       else "Shape of passed values")
                                self.assertRaisesRegex(ValueError,
                                                       err,
                                                       reframe,
                                                       value,
                                                       (nrow, cols))
                        else:
                            self.assertRaisesRegex(ValueError,
                                                   "Must pass 2-d input",
                                                   reframe,
                                                   value,
                                                   (nrow, cols))

    def test_array_index_index(self):
        for ndim in range(4):
            for shape in product(range(4), repeat=ndim):
                self.assertIsInstance(shape, tuple)
                self.assertEqual(len(shape), ndim)
                value = rng.random(shape)
                self.assertIsInstance(value, np.ndarray)
                self.assertEqual(value.shape, shape)
                for nrow in range(5):
                    rows = rng.integers(10, size=nrow)
                    self.assertEqual(rows.shape, (nrow,))
                    for ncol in range(3):
                        cols = rng.integers(10, size=ncol)
                        self.assertEqual(cols.shape, (ncol,))
                        if 1 <= ndim <= 2:
                            if ((ndim, nrow, ncol) == (1, shape[0], 1)
                                    or (nrow, ncol) == shape):
                                frame = reframe(value, (rows, cols))
                                self.assertIsInstance(frame, pd.DataFrame)
                                self.assertEqual(frame.shape, (nrow, ncol))
                                if ndim == 1:
                                    self.assertTrue(np.allclose(frame.T, value))
                                else:
                                    self.assertTrue(np.allclose(frame, value))
                                self.assertIsInstance(frame.index, pd.Index)
                                self.assertTrue(np.all(frame.index == rows))
                                self.assertIsInstance(frame.columns, pd.Index)
                                self.assertTrue(np.all(frame.columns == cols))
                            else:
                                err = ("Empty data passed with indices"
                                       if shape[0] == 0 and nrow > 0
                                       else "Shape of passed values")
                                self.assertRaisesRegex(ValueError,
                                                       err,
                                                       reframe,
                                                       value,
                                                       (rows, cols))
                        else:
                            self.assertRaisesRegex(ValueError,
                                                   "Must pass 2-d input",
                                                   reframe,
                                                   value,
                                                   (rows, cols))


if __name__ == "__main__":
    ut.main()

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
