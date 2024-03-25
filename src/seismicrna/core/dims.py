from typing import Iterable, Sequence

import numpy as np


def find_dims(dims: Sequence[Sequence[str | None]],
              arrays: Sequence[np.ndarray],
              names: Sequence[str] | None = None,
              nonzero: Iterable[str] | bool = False):
    """ Check the dimensions of the arrays.

    Parameters
    ----------

    """
    # Ensure that nonzero is either True or a list of str.
    if nonzero is False:
        nonzero = list()
    elif nonzero is not True:
        nonzero = list(map(str, nonzero))
    # Verify there are the same number of arrays, dimensions, and names.
    if (n := len(arrays)) != len(dims):
        raise ValueError("The numbers of arrays and dimensions must equal, "
                         f"but got {n} array(s) and {len(dims)} dimension(s)")
    if names is not None:
        if len(names) != n:
            raise ValueError("The numbers of arrays and names must equal, "
                             f"but got {n} array(s) and {len(names)} name(s)")
    else:
        names = [f"array{i}" for i in range(n)]
    # Check the dimensions of the arrays.
    sizes = dict()
    for array, dim, name in zip(arrays, dims, names, strict=True):
        if not isinstance(array, np.ndarray):
            raise TypeError(f"Each array must be a NumPy NDArray, "
                            f"but got {type(array).__name__} for {repr(name)}")
        # Count the named and extra dimensions for this array.
        if len(dim) > 0 and dim[-1] is None:
            # The last dimension is None, so subtract it from the number
            # of named dimensions.
            n_named = len(dim) - 1
            # Extra dimensions in the array are allowed.
            extras = True
        else:
            # All dimensions are named.
            n_named = len(dim)
            # Extra dimensions in the array are forbidden.
            extras = False
        # Verify the array has a valid number of dimensions.
        if array.ndim != n_named:
            if not extras:
                raise ValueError(f"Array {repr(name)} must have {n_named} "
                                 f"dimension(s), but got {array.ndim}")
            if array.ndim < n_named:
                raise ValueError(f"Array {repr(name)} must have ≥ {n_named} "
                                 f"dimension(s), but got {array.ndim}")
        # Check each named dimension of the array.
        for i in range(n_named):
            if not isinstance(dim[i], str):
                raise TypeError("The name of each dimension must be str, "
                                f"but got {type(dim[i]).__name__}")
            # Get the size of this dimension in the array.
            size = array.shape[i]
            if (other_size := sizes.get(dim[i])) is not None:
                # A dimension of this name was already encountered.
                if size != other_size:
                    raise ValueError("Got multiple sizes for dimension "
                                     f"{repr(dim[i])}: {other_size} ≠ {size}")
            else:
                # This is the first time this dimension was encountered.
                # Validate the size.
                if not isinstance(size, int):
                    raise TypeError(f"Size of dimension {repr(dim[i])} must "
                                    f"be int, but got {type(size).__name__}")
                if nonzero is True:
                    min_size = 1
                elif nonzero is False:
                    min_size = 0
                else:
                    min_size = 1 if dim[i] in nonzero else 0
                if size < min_size:
                    raise ValueError(f"Size of dimension {repr(dim[i])} must "
                                     f"be ≥ {min_size}, but got {size}")
                sizes[dim[i]] = size
    # Check if any dimensions in nonzero were not defined.
    if nonzero is not True:
        for dim in nonzero:
            if dim not in sizes:
                raise ValueError(f"Unknown dimension for nonzero: {repr(dim)}")
    # Return the size of each dimension.
    return sizes


def triangular(n: int):
    """ The `n` th triangular number (`n` ≥ 0): number of items in an
    equilateral triangle with `n` items on each side.

    Parameters
    ----------
    n: int
        Index of the triangular number to return; equivalently, the side
        length of the equilateral triangle.

    Returns
    -------
    int
        The triangular number with index `n`; equivalently, the number
        of items in the equilateral triangle of side length `n`.
    """
    if not isinstance(n, int):
        raise TypeError(f"n must be int, but got {type(n).__name__}")
    if n < 0:
        raise ValueError(f"n must be ≥ 0, but got {n}")
    return (n * n + n) // 2
