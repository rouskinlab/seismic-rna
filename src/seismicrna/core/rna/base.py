from functools import cached_property

from ..seq import InvalidBaseError, DNA, RNA, Region


class RNARegion(object):
    """ Region of an RNA sequence. """

    def __init__(self, *, region: Region, **kwargs):
        super().__init__(**kwargs)
        self.region = region

    @property
    def init_args(self):
        """ Arguments needed to initialize a new instance. """
        return dict(region=self.region)

    @property
    def ref(self):
        """ Name of the reference sequence. """
        return self.region.ref

    @property
    def end5(self):
        """ Position of the 5' end of the region. """
        return self.region.end5

    @property
    def end3(self):
        """ Position of the 3' end of the region. """
        return self.region.end3

    @property
    def reg(self):
        """ Name of the region. """
        return self.region.name

    @cached_property
    def seq(self):
        """ Sequence of the region as RNA. """
        return self.region.seq.tr()

    @property
    def seq_record(self):
        return self.region.ref_reg, self.seq

    def _subregion_kwargs(self, end5: int, end3: int):
        """ Keyword arguments used by self.subregion(). """
        return dict(region=self.region.subregion(end5, end3))

    def subregion(self, end5: int, end3: int):
        return self.__class__(**self._subregion_kwargs(end5, end3))

    def _renumber_from_args(self, seq5: int):
        """ Arguments needed to initialize a renumbered instance. """
        return self.init_args | dict(region=self.region.renumber_from(seq5))

    def renumber_from(self, seq5: int):
        """ Return a new RNARegion renumbered starting from a position.

        Parameters
        ----------
        seq5: int
            Position from which to start the new numbering system.

        Returns
        -------
        RNARegion
            RNARegion with renumbered positions.
        """
        return self.__class__(**self._renumber_from_args(seq5))

    def __str__(self):
        return f"{type(self).__name__} over {self.region}"
