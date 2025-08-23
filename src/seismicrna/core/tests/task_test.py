import unittest as ut

from seismicrna.core.task import calc_pool_size


class TestCalcPoolSize(ut.TestCase):

    def test_1_task(self):
        for num_cpus in range(1, 5):
            with self.subTest(num_cpus=num_cpus):
                expect = 1, num_cpus
                self.assertTupleEqual(
                    calc_pool_size(1, num_cpus),
                    expect
                )

    def test_1_proc(self):
        for n_tasks in range(1, 5):
            with self.subTest(n_tasks=n_tasks):
                expect = 1, 1
                self.assertTupleEqual(
                    calc_pool_size(n_tasks, 1),
                    expect
                )

    def test_multiple(self):
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             num_cpus=2),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             num_cpus=3),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             num_cpus=4),
                              (2, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             num_cpus=5),
                              (2, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             num_cpus=2),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             num_cpus=3),
                              (3, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             num_cpus=6),
                              (3, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             num_cpus=7),
                              (3, 2))


if __name__ == "__main__":
    ut.main(verbosity=2)
