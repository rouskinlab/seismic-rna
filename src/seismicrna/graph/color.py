from abc import ABC, abstractmethod
from functools import cache
from inspect import getmembers
from sys import modules
from typing import Any, Hashable

from ..core.seq import A_INT, C_INT, G_INT, T_INT
from ..table.base import REL_CODES


class ColorMap(ABC):
    """ Color map for a graph. """

    def __init__(self, name: str, **kwargs):
        self.name = name
        self._colors = self._set_colors(**kwargs)

    @abstractmethod
    def _set_colors(self, **kwargs):
        return dict(**kwargs)

    def get(self, item: Hashable, default: Any | None = None):
        return self._colors.get(item, default)

    def __getitem__(self, item: Hashable):
        return self._colors[item]


class SeqColorMap(ColorMap):
    """ Color map for bases A, C, G, and T. """

    def __init__(self, name: str, a: str, c: str, g: str, t: str):
        super().__init__(name, a=a, c=c, g=g, t=t)

    def _set_colors(self, *, a: str, c: str, g: str, t: str):
        return {A_INT: a, C_INT: c, G_INT: g, T_INT: t}


class RelColorMap(ColorMap):
    """ Color map for relationships. """

    def __init__(self, name: str, v: str, n: str, r: str, m: str,
                 d: str, i: str, s: str, a: str, c: str, g: str, t: str):
        super().__init__(name, v=v, n=n, r=r, m=m,
                         d=d, i=i, s=s, a=a, c=c, g=g, t=t)

    def _set_colors(self, **kwargs):
        colors = dict()
        for key, field in REL_CODES.items():
            colors[field] = kwargs.pop(key)
        if kwargs:
            raise TypeError(f"Unexpected keyword arguments: {kwargs}")
        return colors


basic = SeqColorMap("basic", a="#FF0000", c="#0000FF", g="#FFC000", t="#008000")
water = SeqColorMap("water", a="#A15252", c="#3D427D", g="#E3CC7B", t="#76B887")
earth = SeqColorMap("earth", a="#D17777", c="#464EA6", g="#E3CC7B", t="#336140")
steel = SeqColorMap("steel", a="#663328", c="#716B80", g="#91B8AC", t="#D9D5B4")
tetra = SeqColorMap("tetra", a="#C05F15", c="#597DE4", g="#743B4A", t="#9BD1D0")

crayons = RelColorMap("crayons", v="#424242", n="#A9A9A9", r="#942193",
                      m="#929000", d="#FF2600", i="#00FA92", s="#FF40FF",
                      a="#73FCD6", c="#FFD479", g="#7A81FF", t="#FF8AD8")
sexta = RelColorMap("sexta", v="#FBED94", n="#C05F15", r="#597DE4",
                    m="#0F155F", d="#FBED94", i="#0F155F", s="#743B4A",
                    a="#C05F15", c="#597DE4", g="#743B4A", t="#9BD1D0")

DEFAULTS: dict[type[ColorMap], ColorMap] = {
    RelColorMap: sexta,
    SeqColorMap: tetra,
}


@cache
def get_colormaps(cmap_class: type[ColorMap]):
    """ Return a dict of all color maps of a given class. """
    colormaps: dict[str, cmap_class] = dict()
    for _, cmap in getmembers(modules[__name__],
                              lambda item: isinstance(item, cmap_class)):
        if cmap.name in colormaps:
            raise ValueError(f"Duplicate {cmap_class.__name__}: '{cmap.name}'")
        colormaps[cmap.name] = cmap
    if (default := DEFAULTS[cmap_class]) not in colormaps:
        raise ValueError(
            f"Default {cmap_class.__name__} '{default}' does not exist")
    return colormaps


def get_cmap(cmap_class: type[ColorMap], name: str | None = None):
    """ Get a color map of a given class by its name. """
    if name is None:
        # Use the default color map for the class.
        return DEFAULTS[cmap_class]
    cmaps = get_colormaps(cmap_class)
    if not cmaps:
        raise ValueError(f"No color maps of class {cmap_class.__name__}")
    return cmaps[name]
