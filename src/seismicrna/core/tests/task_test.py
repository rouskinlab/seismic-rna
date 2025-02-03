import unittest as ut

from seismicrna.core.task import calc_pool_size


class TestCalcPoolSize(ut.TestCase):

    def test_1_task(self):
        for max_procs in range(1, 5):
            with self.subTest(max_procs=max_procs):
                expect = 1, max_procs
                self.assertTupleEqual(
                    calc_pool_size(1, max_procs),
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
                                             max_procs=2),
                              (1, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=3),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=4),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=5),
                              (2, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=2),
                              (1, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=3),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=6),
                              (3, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=7),
                              (3, 2))


if __name__ == "__main__":
    ut.main(verbosity=2)
