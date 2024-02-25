from ..core import path
from ..table.base import Table


def get_list_path(table: Table):
    """ Get the path to a list file corresponding to a table. """
    path_fields = table.path_fields
    path_fields[path.CMD] = path.CMD_LIST_DIR
    return path.buildpar(*path.POS_TABLE_SEGS, **path_fields)
