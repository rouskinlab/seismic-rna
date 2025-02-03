import unittest as ut

import numpy as np

from seismicrna.sim.ends import sim_pends


class TestSimPEnds(ut.TestCase):

    def equate(self,
               end5: int,
               end3: int,
               center_fmean: float,
               center_fvar: float,
               length_fmean: float,
               length_fvar: float,
               end5s: np.ndarray,
               end3s: np.ndarray,
               pends: np.ndarray,
               **kwargs):
        for result, expect in zip(sim_pends(end5,
                                            end3,
                                            center_fmean,
                                            center_fvar,
                                            length_fmean,
                                            length_fvar,
                                            **kwargs),
                                  [end5s, end3s, pends],
                                  strict=True):
            self.assertTrue(np.array_equal(result, expect))

    def test_zero_length_1_to_0_keep(self):
        self.equate(end5=1,
                    end3=0,
                    center_fmean=0.5,
                    center_fvar=0.5,
                    length_fmean=0.5,
                    length_fvar=0.5,
                    end5s=np.array([1]),
                    end3s=np.array([0]),
                    pends=np.array([1.]))

    def test_zero_length_9_to_8_keep(self):
        self.equate(end5=9,
                    end3=8,
                    center_fmean=0.5,
                    center_fvar=0.5,
                    length_fmean=0.5,
                    length_fvar=0.5,
                    end5s=np.array([9]),
                    end3s=np.array([8]),
                    pends=np.array([1.]))

    def test_zero_length_9_to_8_drop(self):
        self.equate(end5=9,
                    end3=8,
                    center_fmean=0.5,
                    center_fvar=0.5,
                    length_fmean=0.5,
                    length_fvar=0.5,
                    end5s=np.array([]),
                    end3s=np.array([]),
                    pends=np.array([]),
                    keep_empty_reads=False)

    def test_center_fvar_0_length_fvar_0(self):
        self.equate(end5=11,
                    end3=19,
                    center_fmean=0.,
                    center_fvar=0.,
                    length_fmean=1.,
                    length_fvar=0.,
                    end5s=np.array([11]),
                    end3s=np.array([11]),
                    pends=np.array([1.]))
        self.equate(end5=11,
                    end3=19,
                    center_fmean=0.25,
                    center_fvar=0.,
                    length_fmean=1.,
                    length_fvar=0.,
                    end5s=np.array([11]),
                    end3s=np.array([15]),
                    pends=np.array([1.]))
        self.equate(end5=11,
                    end3=19,
                    center_fmean=0.5,
                    center_fvar=0.,
                    length_fmean=1.,
                    length_fvar=0.,
                    end5s=np.array([11]),
                    end3s=np.array([19]),
                    pends=np.array([1.]))
        self.equate(end5=11,
                    end3=19,
                    center_fmean=0.5,
                    center_fvar=0.,
                    length_fmean=0.5,
                    length_fvar=0.,
                    end5s=np.array([13]),
                    end3s=np.array([17]),
                    pends=np.array([1.]))
        self.equate(end5=11,
                    end3=19,
                    center_fmean=0.75,
                    center_fvar=0.,
                    length_fmean=1.,
                    length_fvar=0.,
                    end5s=np.array([15]),
                    end3s=np.array([19]),
                    pends=np.array([1.]))
        self.equate(end5=11,
                    end3=19,
                    center_fmean=1.,
                    center_fvar=0.,
                    length_fmean=1.,
                    length_fvar=0.,
                    end5s=np.array([19]),
                    end3s=np.array([19]),
                    pends=np.array([1.]))


if __name__ == "__main__":
    ut.main(verbosity=2)
