from ..core import path
from ..core.table import Table


def get_list_path(table: Table):
    """ Get the path to a list file corresponding to a table. """
    path_fields = table.get_default_path_fields() | table.get_path_fields()
    path_fields[path.CMD] = path.LIST_STEP
    return path.buildpar(*table.get_path_segs(), **path_fields)
