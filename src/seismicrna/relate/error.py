class RelateError(Exception):
    """ Any error that occurs during relating. """


class RelateValueError(RelateError, ValueError):
    """ Any ValueError that occurs during relating. """


class RelateNotImplementedError(RelateError, NotImplementedError):
    """ Any NotImplementedError that occurs during relating. """
