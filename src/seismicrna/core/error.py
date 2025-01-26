""" Generic Exceptions """


class IncompatibleValuesError(ValueError):
    """ Two or more values are individually valid, but their combination
    is not. """


class InconsistentValueError(ValueError):
    """ Two or more values differ when they should be equal. """


class OutOfBoundsError(ValueError):
    """ A numeric value is outside its proper bounds. """
