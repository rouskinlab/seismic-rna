from itertools import filterfalse, tee
from typing import Any, Callable, Iterable


def filterboth(function: Callable[[Any], bool] | None, iterable: Iterable):
    """ Return an iterator yielding those items of iterable for which
    function(item) is true. If function is None, return the items that
    are true. This function is based on a recipe from itertools:
    https://docs.python.org/3/library/itertools.html#itertools-recipes

    Parameters
    ----------
    function: Callable[[Any], bool] | None
        Function to call with each element to determine its truthiness.
        If None, then use the truthiness of the element itself.
    iterable: Iterable
        Items to filter.

    Returns
    -------
    tuple[Iterable, Iterable]
        Two iterables containing the true and false items, respectively.
    """
    i1, i2 = tee(iterable, 2)
    return filter(function, i2), filterfalse(function, i1)


def partition(function: Callable[[Any], bool] | None, iterable: Iterable):
    """ Similar to filterboth, but returns two lists, not iterators. """
    t, f = filterboth(function, list(iterable))
    return list(t), list(f)
