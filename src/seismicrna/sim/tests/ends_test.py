import unittest as ut

from seismicrna.sim.ends import _sim_ends


class TestSimEnds(ut.TestCase):

    def test_typical(self):
        num_reads = 100_000
        for end5 in [1, 11]:
            for end3 in [90, 100]:
                for mean3 in [60., 80.]:
                    for meanr in [20., 40.]:
                        for var in [0.1, 0.5]:
                            kwargs = dict(end5=end5,
                                          end3=end3,
                                          end3_mean=mean3,
                                          read_mean=meanr,
                                          variance=var)
                            with self.subTest(**kwargs):
                                end5s, end3s = _sim_ends(num_reads=num_reads,
                                                         **kwargs)
                                lengths = end3s - end5s + 1
                                mean5 = mean3 - meanr + 1
                                self.assertEqual(round(end5s.mean()), mean5)
                                self.assertEqual(round(end3s.mean()), mean3)
                                self.assertEqual(round(lengths.mean()), meanr)

    def test_no_gap5(self):
        num_reads = 100_000
        for end5 in [1, 11]:
            for end3 in [90, 100]:
                for mean3 in [60., 80.]:
                    meanr = mean3 - (end5 - 1)
                    for var in [0.1, 0.5]:
                        kwargs = dict(end5=end5,
                                      end3=end3,
                                      end3_mean=mean3,
                                      read_mean=meanr,
                                      variance=var)
                        with self.subTest(**kwargs):
                            end5s, end3s = _sim_ends(num_reads=num_reads,
                                                     **kwargs)
                            lengths = end3s - end5s + 1
                            mean5 = mean3 - meanr + 1
                            self.assertEqual(round(end5s.mean()), mean5)
                            self.assertEqual(round(end3s.mean()), mean3)
                            self.assertEqual(round(lengths.mean()), meanr)

    def test_read_zero(self):
        num_reads = 100_000
        meanr = 0.
        for end5 in [1, 11]:
            for end3 in [90, 100]:
                for mean3 in [60., 80.]:
                    for var in [0.1, 0.5]:
                        kwargs = dict(end5=end5,
                                      end3=end3,
                                      end3_mean=mean3,
                                      read_mean=meanr,
                                      variance=var)
                        with self.subTest(**kwargs):
                            end5s, end3s = _sim_ends(num_reads=num_reads,
                                                     **kwargs)
                            lengths = end3s - end5s + 1
                            mean5 = mean3 - meanr + 1
                            self.assertEqual(round(end5s.mean()), mean5)
                            self.assertEqual(round(end3s.mean()), mean3)
                            self.assertEqual(round(lengths.mean()), meanr)

    def test_no_gap3(self):
        num_reads = 100_000
        for end5 in [1, 11]:
            for end3 in [90, 100]:
                mean3 = end3
                for meanr in [20., 40.]:
                    for var in [0.1, 0.5]:
                        kwargs = dict(end5=end5,
                                      end3=end3,
                                      end3_mean=mean3,
                                      read_mean=meanr,
                                      variance=var)
                        with self.subTest(**kwargs):
                            end5s, end3s = _sim_ends(num_reads=num_reads,
                                                     **kwargs)
                            lengths = end3s - end5s + 1
                            mean5 = mean3 - meanr + 1
                            self.assertEqual(round(end5s.mean()), mean5)
                            self.assertEqual(round(end3s.mean()), mean3)
                            self.assertEqual(round(lengths.mean()), meanr)

    def test_var_zero(self):
        num_reads = 100_000
        var = 0.
        for end5 in [1, 11]:
            for end3 in [90, 100]:
                for mean3 in [60., 80.]:
                    for meanr in [20., 40.]:
                        kwargs = dict(end5=end5,
                                      end3=end3,
                                      end3_mean=mean3,
                                      read_mean=meanr,
                                      variance=var)
                        with self.subTest(**kwargs):
                            end5s, end3s = _sim_ends(num_reads=num_reads,
                                                     **kwargs)
                            lengths = end3s - end5s + 1
                            mean5 = mean3 - meanr + 1
                            self.assertEqual(round(end5s.mean()), mean5)
                            self.assertEqual(round(end3s.mean()), mean3)
                            self.assertEqual(round(lengths.mean()), meanr)


if __name__ == "__main__":
    ut.main()
