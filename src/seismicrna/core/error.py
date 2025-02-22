""" Generic Exceptions """


class IncompatibleValuesError(ValueError):
    """ Two or more values are individually valid, but their combination
    is not. """


class IncompatibleOptionsError(ValueError):
    """ Two or more options are incompatible. """


class InconsistentValueError(ValueError):
    """ Two or more values differ when they should be equal. """


class DuplicateValueError(ValueError):
    """ A value occurred more than once when all must be unique. """


class OutOfBoundsError(ValueError):
    """ A numeric value is outside its proper bounds. """


class NoDataError(RuntimeError):
    """ Data were required, but none were provided. """
