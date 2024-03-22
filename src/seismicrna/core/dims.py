from typing import Sequence

import numpy as np


def find_dims(dims: Sequence[Sequence[str | None]],
              arrays: Sequence[np.ndarray],
              names: Sequence[str] | None = None):
    """ Check the dimensions of the arrays.

    Parameters
    ----------

    """
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
                sizes[dim[i]] = size
    # Return the size of each dimension.
    return sizes
