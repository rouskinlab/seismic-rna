from inspect import getmembers
from sys import modules
from typing import Iterable


def get_subclasses(types: type | tuple[type], items: Iterable):
    """ Yield every item from `items` that is a subclass of `types`. """
    return [i for i in items if isinstance(i, type) and issubclass(i, types)]


def get_subclasses_members(types: type | tuple[type], item: object):
    """ Yield every member of `item` that is a subclass of `types`. """
    return get_subclasses(types, (value for name, value in getmembers(item)))


def get_subclasses_module(types: type | tuple[type], module_name: str):
    """ Yield every member of the module named `module_name` that is a
    subclass of `types`. """
    return get_subclasses_members(types, modules[module_name])
