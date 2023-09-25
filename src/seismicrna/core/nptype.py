import numpy as np

# Constants
BITS_PER_BYTE = 8
UINT_CODE = 'u'
BYTE_CODE = 'S'
ENDIAN_STR = "little"
ENDIAN_TYP = '<'


def get_dtype(code: str, size: int):
    """ NumPy type with the given code and size. """
    return np.dtype(f"{ENDIAN_TYP}{code}{size}")


def get_uint_dtype(nbytes: int):
    """ NumPy uint data type with the given number of bytes. """
    return get_dtype(UINT_CODE, nbytes)


def get_uint_type(nbytes: int):
    """ NumPy uint type with the given number of bytes. """
    return get_uint_dtype(nbytes).type


def get_uint_size(uint_type: type):
    """ Size of a NumPy uint type in bytes. """
    return uint_type(0).itemsize


def get_uint_max(uint_type: type):
    """ Maximum value of a NumPy data type. """
    return 2 ** (BITS_PER_BYTE * get_uint_size(uint_type)) - 1


def get_byte_dtype(nchars: int):
    """ NumPy byte type with the given number of characters. """
    return get_dtype(BYTE_CODE, nchars)
