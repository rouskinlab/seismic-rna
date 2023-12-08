import numpy as np
import pandas as pd


def recontainerize(values: np.ndarray | pd.Series | pd.DataFrame,
                   container: np.ndarray | pd.Series | pd.DataFrame):
    """ Return the given values in a container of the given type and,
    if applicable, with the given indexes and columns.

    Parameters
    ----------
    values: numpy.ndarray | pandas.Series | pandas.DataFrame
        Values to copy into the new container.
    container: numpy.ndarray | pandas.Series | pandas.DataFrame
        Type of the container to return; indexes and columns will be
        returned in the output as well.

    Returns
    -------
    numpy.ndarray | pandas.Series | pandas.DataFrame
        Values in the new container.
    """
    if values.shape != container.shape:
        raise ValueError(f"Shapes differ between values {values.shape} and "
                         f"container {container.shape}")
    if isinstance(container, pd.DataFrame):
        return pd.DataFrame(values, container.index, container.columns)
    if isinstance(container, pd.Series):
        return pd.Series(values, container.index)
    if isinstance(container, np.ndarray):
        return np.asarray(values)
    raise TypeError(f"Invalid type for container: {type(container).__name__}")
