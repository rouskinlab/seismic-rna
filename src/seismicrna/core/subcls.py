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
