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
    # Compute the difference in the dimensionality of values and target.
    ddim = len(tshape) - len(vshape)
    if ddim < 0:
        # An array cannot be broadcast to a shape with fewer dimensions.
        return False
    if ddim > 0:
        # If the target shape has more dimensions than the array, then
        # prepend extra dimensions (each of length 1) to the array.
        return broadcastable(ddim * (1,) + vshape, tshape)
    # The array is broadcastable to the target shape if and only if, for
    # every dimension, the length of the original axis is 1 or equals
    # the length of the target axis.
    for vlen, tlen in zip(vshape, tshape, strict=True):
        if vlen != 1 and vlen != tlen:
            return False
    return True


def broadcast_regex(vshape: tuple[int, ...], tshape: tuple[int, ...]):
    """ Get the error regex that would result from trying to broadcast
    values in the shape of `vshape` to the target shape of `tshape`.

    Parameters
    ----------
    vshape: tuple[int, ...]
        Shape of the values to broadcast.
    tshape: tuple[int, ...]
        Shape of the target to which the values should be broadcast.

    Returns
    -------
    str
        Regular expression that can match the error message.
    """
    if broadcastable(vshape, tshape):
        raise ValueError(f"Broadcasting {vshape} to {tshape} is allowed")
    vdim = len(vshape)
    tdim = len(tshape)
    if tdim == 0:
        return "cannot broadcast a non-scalar to a scalar array"
    elif vdim > tdim:
        return "input operand has more dimensions than allowed"
    else:
        return "could not be broadcast together"


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
            for ncol in range(3):
                cols = rng.integers(10, size=ncol)
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
            for ncol in range(3):
                cols = rng.integers(10, size=ncol)
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
            for ncol in range(3):
                cols = rng.integers(10, size=ncol)
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
                value = rng.random(shape)
                frame = reframe(value, shape)
                self.assertIsInstance(frame, np.ndarray)
                self.assertEqual(frame.shape, shape)
                self.assertTrue(np.allclose(frame, value))

    def test_array_ints(self):
        for vdim in range(4):
            for vshape in product(range(4), repeat=vdim):
                value = rng.random(vshape)
                for tdim in range(4):
                    for tshape in product(range(4), repeat=tdim):
                        if broadcastable(vshape, tshape):
                            fshape = np.broadcast_shapes(vshape, tshape)
                            frame = reframe(value, tshape)
                            self.assertIsInstance(frame, np.ndarray)
                            self.assertEqual(frame.shape, fshape)
                            self.assertTrue(np.allclose(frame, value))
                        else:
                            self.assertRaisesRegex(ValueError,
                                                   broadcast_regex(vshape,
                                                                   tshape),
                                                   reframe,
                                                   value,
                                                   tshape)

    def test_array_index(self):
        for ndim in range(4):
            for shape in product(range(4), repeat=ndim):
                value = rng.random(shape)
                for length in range(5):
                    index = rng.integers(10, size=length)
                    if broadcastable(shape, (length,)):
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
                                               broadcast_regex(shape,
                                                               (length,)),
                                               reframe,
                                               value,
                                               (index,))

    def test_array_index_index(self):
        for ndim in range(4):
            for vshape in product(range(4), repeat=ndim):
                value = rng.random(vshape)
                for nrow in range(5):
                    rows = rng.integers(10, size=nrow)
                    for ncol in range(3):
                        tshape = nrow, ncol
                        cols = rng.integers(10, size=ncol)
                        axes_sets = [(rows, ncol), (rows, cols), (nrow, cols)]
                        for axes in axes_sets:
                            if broadcastable(vshape, tshape):
                                frame = reframe(value, axes)
                                self.assertIsInstance(frame, pd.DataFrame)
                                self.assertEqual(frame.shape, tshape)
                                self.assertTrue(np.allclose(frame, value))
                                if isinstance(axes[0], int):
                                    self.assertIsInstance(frame.index,
                                                          pd.RangeIndex)
                                    self.assertTrue(np.all(frame.index
                                                           == np.arange(nrow)))
                                else:
                                    self.assertIsInstance(frame.index,
                                                          pd.Index)
                                    self.assertTrue(np.all(frame.index
                                                           == rows))
                                if isinstance(axes[1], int):
                                    self.assertIsInstance(frame.columns,
                                                          pd.RangeIndex)
                                    self.assertTrue(np.all(frame.columns
                                                           == np.arange(ncol)))
                                else:
                                    self.assertIsInstance(frame.columns,
                                                          pd.Index)
                                    self.assertTrue(np.all(frame.columns
                                                           == cols))
                            else:
                                self.assertRaisesRegex(ValueError,
                                                       broadcast_regex(vshape,
                                                                       tshape),
                                                       reframe,
                                                       value,
                                                       axes)

    def test_series(self):
        for vlength in range(5):
            vshape = vlength,
            value = pd.Series(rng.random(vlength))
            for tlength in range(5):
                tshape = tlength,
                index = rng.integers(10, size=tlength)
                axes_sets = [(index,), tshape, None]
                for axes in axes_sets:
                    if axes is None or broadcastable(vshape, tshape):
                        frame = reframe(value, axes)
                        if axes is None:
                            self.assertIsInstance(frame, np.ndarray)
                            self.assertEqual(frame.shape, vshape)
                            self.assertTrue(np.allclose(frame, value))
                        elif all(isinstance(axis, int) for axis in axes):
                            self.assertIsInstance(frame, np.ndarray)
                            self.assertEqual(frame.shape, tshape)
                            self.assertTrue(np.allclose(frame, value))
                        else:
                            self.assertIsInstance(frame, pd.Series)
                            self.assertEqual(frame.shape, tshape)
                            self.assertTrue(np.allclose(frame, value))
                            self.assertIsInstance(frame.index, pd.Index)
                            self.assertTrue(np.all(frame.index == index))
                    else:
                        self.assertRaisesRegex(ValueError,
                                               broadcast_regex(vshape,
                                                               tshape),
                                               reframe,
                                               value,
                                               axes)

    def test_dataframe(self):
        for vshape in product(range(4), repeat=2):
            value = pd.DataFrame(rng.random(vshape))
            for nrow in range(5):
                rows = rng.integers(10, size=nrow)
                for ncol in range(3):
                    tshape = nrow, ncol
                    cols = rng.integers(10, size=ncol)
                    axes_sets = [(rows, cols),
                                 (rows, ncol),
                                 (nrow, cols),
                                 (nrow, ncol),
                                 None]
                    for axes in axes_sets:
                        if axes is None or broadcastable(vshape, tshape):
                            frame = reframe(value, axes)
                            if axes is None:
                                self.assertIsInstance(frame, np.ndarray)
                                self.assertEqual(frame.shape, vshape)
                                self.assertTrue(np.allclose(frame, value))
                            elif all(isinstance(axis, int) for axis in axes):
                                self.assertIsInstance(frame, np.ndarray)
                                self.assertEqual(frame.shape, tshape)
                                self.assertTrue(np.allclose(frame, value))
                            else:
                                self.assertIsInstance(frame, pd.DataFrame)
                                self.assertEqual(frame.shape, tshape)
                                self.assertTrue(np.allclose(frame, value))
                                if isinstance(axes[0], int):
                                    self.assertIsInstance(frame.index,
                                                          pd.RangeIndex)
                                    self.assertTrue(np.all(frame.index
                                                           == np.arange(nrow)))
                                else:
                                    self.assertIsInstance(frame.index,
                                                          pd.Index)
                                    self.assertTrue(np.all(frame.index
                                                           == rows))
                                if isinstance(axes[1], int):
                                    self.assertIsInstance(frame.columns,
                                                          pd.RangeIndex)
                                    self.assertTrue(np.all(frame.columns
                                                           == np.arange(ncol)))
                                else:
                                    self.assertIsInstance(frame.columns,
                                                          pd.Index)
                                    self.assertTrue(np.all(frame.columns
                                                           == cols))
                        else:
                            self.assertRaisesRegex(ValueError,
                                                   broadcast_regex(vshape,
                                                                   tshape),
                                                   reframe,
                                                   value,
                                                   axes)


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
