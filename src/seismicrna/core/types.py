from functools import cache

import numpy as np

# Constants
BITS_PER_BYTE = 8
UINT_CODE = "u"
BYTE_CODE = "S"
ENDIAN_TYP = "<"
UINT_NBYTES = 1, 2, 4, 8


def get_dtype(code: str, size: int):
    """ NumPy type with the given code and size. """
    return np.dtype(f"{ENDIAN_TYP}{code}{size}")


@cache
def get_uint_dtype(nbytes: int):
    """ NumPy uint data type with the given number of bytes. """
    if nbytes not in UINT_NBYTES:
        raise ValueError(f"Invalid number of bytes for unsigned int: {nbytes}")
    return get_dtype(UINT_CODE, nbytes)


def get_uint_type(nbytes: int):
    """ NumPy uint type with the given number of bytes. """
    return get_uint_dtype(nbytes).type


@cache
def get_uint_size(uint_type: type):
    """ Size of a NumPy uint type in bytes. """
    return uint_type(0).itemsize


@cache
def get_max_value(nbytes: int):
    """ Get the maximum value of an unsigned integer of N bytes. """
    return 2 ** (BITS_PER_BYTE * nbytes) - 1


def get_max_uint(uint_type: type):
    """ Maximum value of a NumPy unsigned integer type. """
    return get_max_value(get_uint_size(uint_type))


@cache
def fit_uint_size(value: int):
    """ Smallest number of bytes that will fit the value. """
    if not isinstance(value, (int, np.integer)) or value < 0:
        raise ValueError(f"Expected an integer ≥ 0, but got {repr(value)}")
    for nbytes in UINT_NBYTES:
        if value <= get_max_value(nbytes):
            return nbytes
    raise ValueError(f"Value does not fit into any integer type: {value}")


def fit_uint_type(value: int):
    """ Smallest unsigned int type that will fit the value. """
    return get_uint_type(fit_uint_size(value))


@cache
def get_byte_dtype(nchars: int):
    """ NumPy byte type with the given number of characters. """
    return get_dtype(BYTE_CODE, nchars)

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
