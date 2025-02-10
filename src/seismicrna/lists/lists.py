from ..core import path
from ..core.table import Table


def get_list_path(table: Table):
    """ Get the path to a list file corresponding to a table. """
    path_fields = table.get_path_field_values()
    path_fields[path.STEP] = path.LIST_STEP
    return path.buildpar(*table.get_path_seg_types(), **path_fields)
