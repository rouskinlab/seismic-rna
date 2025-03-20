from operator import eq, ne, ge, gt, le, lt
from typing import Any, Callable, Container, Type

import numpy as np
import pandas as pd


def require_issubclass(name: str,
                       value: type,
                       classes: type | tuple[type | tuple[Any, ...], ...],
                       error_type: Type[ValueError] = ValueError) -> None:
    """ Raise an error if value is not a subclass of classes. """
    if not issubclass(value, classes):
        require_issubclass("error_type", error_type, ValueError)
        raise error_type(f"{name} must be a subclass of {classes}, "
                         f"but got {repr(value)}")


def require_isinstance(name: str,
                       value: Any,
                       classes: type | tuple[type | tuple[Any, ...], ...],
                       error_type: Type[TypeError] = TypeError) -> None:
    """ Raise an error if value is not an instance of classes. """
    if not isinstance(value, classes):
        require_issubclass("error_type", error_type, TypeError)
        raise error_type(f"{name} must be an instance of {classes}, "
                         f"but got {repr(value)} of type {type(value)}")


def require_isin(name: str,
                 value: Any,
                 values: Container,
                 values_name: str = "",
                 error_type: Type[Exception] = ValueError):
    """ Require value to be in values. """
    if value not in values:
        require_issubclass("error_type", error_type, Exception)
        if values_name:
            message = (f"{name} must be in {values_name}, "
                       f"but got {name}={repr(value)} "
                       f"and {values_name}={repr(values)}")
        else:
            message = (f"{name} must be in {repr(values)}, "
                       f"but got {repr(value)}")
        raise error_type(message)


def _require_compare(
        name: str,
        value: Any,
        comparison: Callable[[Any, Any], bool],
        compare_name: str,
        default_compare_name: str,
        compare_value: Any,
        classes: type | tuple[type | tuple[Any, ...], ...],
        error_type: Type[ValueError]
):
    require_isinstance(name, value, classes)
    require_isinstance(compare_name if compare_name else default_compare_name,
                       compare_value,
                       classes)
    # Ensure the comparison function is valid.
    comparison_signs = {eq: "=",
                        ne: "≠",
                        ge: "≥",
                        gt: ">",
                        le: "≤",
                        lt: "<",
                        np.array_equal: "=",
                        np.allclose: "≈",
                        pd.Index.equals: "="}
    require_isin("comparison",
                 comparison,
                 comparison_signs,
                 error_type=(ValueError
                             if callable(comparison)
                             else TypeError))
    if not comparison(value, compare_value):
        require_issubclass("error_type", error_type, ValueError)
        sign = comparison_signs[comparison]
        if compare_name:
            message = (f"Must have {name} {sign} {compare_name}, "
                       f"but got {name}={repr(value)} "
                       f"and {compare_name}={repr(compare_value)}")
        else:
            message = (f"Must have {name} {sign} {repr(compare_value)}, "
                       f"but got {repr(value)}")
        raise error_type(message)


def require_equal(
        name: str,
        value: Any,
        other_value: Any,
        other_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that value = other_value. """
    _require_compare(name, value, eq,
                     other_name, "other", other_value,
                     classes, error_type)


def require_array_equal(
        name: str,
        array: Any,
        other_array: Any,
        other_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that array = other_array. """
    _require_compare(name, array, np.array_equal,
                     other_name, "other", other_array,
                     classes, error_type)


def require_allclose(
        name: str,
        array: Any,
        other_array: Any,
        other_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that array ≈ other_array. """
    _require_compare(name, array, np.allclose,
                     other_name, "other", other_array,
                     classes, error_type)


def require_index_equals(
        name: str,
        index: pd.Index,
        other_index: pd.Index,
        other_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = pd.Index,
        error_type: Type[ValueError] = ValueError
):
    """ Require that index = other_index. """
    _require_compare(name, index, pd.Index.equals,
                     other_name, "other", other_index,
                     classes, error_type)


def require_atleast(
        name: str,
        value: Any,
        minimum_value: Any,
        minimum_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that value ≥ minimum_value. """
    _require_compare(name, value, ge,
                     minimum_name, "minimum", minimum_value,
                     classes, error_type)


def require_greater(
        name: str,
        value: Any,
        other_value: Any,
        other_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that value > other_value. """
    _require_compare(name, value, gt,
                     other_name, "other", other_value,
                     classes, error_type)


def require_atmost(
        name: str,
        value: Any,
        maximum_value: Any,
        maximum_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that value ≤ maximum_value. """
    _require_compare(name, value, le,
                     maximum_name, "maximum", maximum_value,
                     classes, error_type)


def require_less(
        name: str,
        value: Any,
        other_value: Any,
        other_name: str = "",
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that value < other_value. """
    _require_compare(name, value, lt,
                     other_name, "less", other_value,
                     classes, error_type)


def require_between(
        name: str,
        value: Any,
        minimum_value: Any | None,
        maximum_value: Any | None,
        minimum_name: str = "",
        maximum_name: str = "",
        inclusive: bool = True,
        classes: type | tuple[type | tuple[Any, ...], ...] = object,
        error_type: Type[ValueError] = ValueError
):
    """ Require that value is in [minimum_value, maximum_value] if
    inclusive is True, otherwise in (minimum_value, maximum_value). """
    if minimum_value is not None:
        min_args = name, value, minimum_value, minimum_name, classes, error_type
        if inclusive:
            require_atleast(*min_args)
        else:
            require_greater(*min_args)
    if maximum_value is not None:
        max_args = name, value, maximum_value, maximum_name, classes, error_type
        if inclusive:
            require_atmost(*max_args)
        else:
            require_less(*max_args)


def require_fraction(
        name: str,
        value: Any,
        classes: type | tuple[type | tuple[Any, ...], ...] = (float, int),
        error_type: Type[ValueError] = ValueError
):
    """ Require that value ≥ 0 and ≤ 1. """
    require_between(name, value, 0., 1., classes=classes, error_type=error_type)
