from abc import ABC, abstractmethod
from functools import cached_property

from .base import BaseGraph, BaseRunner
from ..core.arg import opt_rels
from ..core.table import get_rel_name


class RelGraph(BaseGraph, ABC):
    """ Graph of one or more types of relationships. """

    @property
    @abstractmethod
    def codestring(self) -> str:
        """ String of the relationship code(s). """

    @property
    @abstractmethod
    def rel_names(self):
        """ Names of the relationships to graph. """

    @cached_property
    def relationships(self) -> str:
        """ Relationships being graphed as a slash-separated string. """
        return "/".join(self.rel_names)


class OneRelGraph(RelGraph, ABC):
    """ Graph of exactly one type of relationship. """

    def __init__(self, *, rel: str, **kwargs):
        """
        Parameters
        ----------
        rel: str
            Relationship(s) whose data will be pulled from the table(s).
            It is given as a one-letter code, the definitions of which
            are given in seismicrna.table.base.REL_CODES, as follows:

            - Covered: v
            - Informative: n
            - Matched: r
            - Mutated: m
            - Subbed: s
            - Subbed to A: a
            - Subbed to C: c
            - Subbed to G: g
            - Subbed to T: t
            - Deleted: d
            - Inserted: i

            This type of graph requires that `rel` be a one-letter code.
        """
        super().__init__(**kwargs)
        if len(rel) != 1:
            raise ValueError(f"{type(self).__name__} expected exactly one "
                             f"relationship code, but got {repr(rel)}")
        self._rel = rel

    @property
    def codestring(self):
        return self._rel

    @cached_property
    def rel_name(self):
        """ Name of the relationship to graph. """
        return get_rel_name(self.codestring)

    @cached_property
    def rel_names(self):
        return [self.rel_name]


class MultiRelsGraph(RelGraph, ABC):
    """ Graph of one or more relationships. """

    def __init__(self, *, rels: str, **kwargs):
        """
        Parameters
        ----------
        rels: str
            Relationships whose data will be pulled from the table(s).
            Each is given as a one-letter code, the definitions of which
            are given in seismicrna.table.base.REL_CODES, as follows:

            - Covered: v
            - Informative: n
            - Matched: r
            - Mutated: m
            - Subbed: s
            - Subbed to A: a
            - Subbed to C: c
            - Subbed to G: g
            - Subbed to T: t
            - Deleted: d
            - Inserted: i

            More than one relationship can be specified by passing a
            string longer than one character; for example, `acgt` would
            mean to use all four types of substitutions.
        """
        super().__init__(**kwargs)
        if len(rels) == 0:
            raise ValueError(f"{type(self).__name__} expected one or more "
                             f"relationship codes, but got {repr(rels)}")
        self._rels = rels

    @property
    def codestring(self):
        return self._rels

    @cached_property
    def rel_names(self):
        return list(map(get_rel_name, self.codestring))


class RelRunner(BaseRunner, ABC):

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_rels]
