from ..core import path
from ..core.table import Table


def get_list_path(table: Table):
    """ Get the path to a list file corresponding to a table. """
    path_fields = table.default_path_fields() | table.path_fields
    path_fields[path.CMD] = path.LIST_STEP
    return path.buildpar(*table.path_segs(), **path_fields)
