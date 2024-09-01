import unittest as ut

from seismicrna.core.task import calc_pool_size


class TestCalcPoolSize(ut.TestCase):

    def test_1_task(self):
        for max_procs in range(1, 5):
            for parallel in [False, True]:
                with self.subTest(max_procs=max_procs, parallel=parallel):
                    expect = 1, max_procs
                    self.assertTupleEqual(
                        calc_pool_size(1, max_procs, parallel),
                        expect
                    )

    def test_1_proc(self):
        for n_tasks in range(1, 5):
            for parallel in [False, True]:
                with self.subTest(n_tasks=n_tasks, parallel=parallel):
                    expect = 1, 1
                    self.assertTupleEqual(
                        calc_pool_size(n_tasks, 1, parallel),
                        expect
                    )

    def test_parallel(self):
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=2,
                                             parallel=True),
                              (1, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=2,
                                             parallel=False),
                              (1, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=3,
                                             parallel=True),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=3,
                                             parallel=False),
                              (1, 3))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=4,
                                             parallel=True),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=4,
                                             parallel=False),
                              (1, 4))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=5,
                                             parallel=True),
                              (2, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=2,
                                             max_procs=5,
                                             parallel=False),
                              (1, 5))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=2,
                                             parallel=True),
                              (1, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=2,
                                             parallel=False),
                              (1, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=3,
                                             parallel=True),
                              (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=3,
                                             parallel=False),
                              (1, 3))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=6,
                                             parallel=True),
                              (3, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=6,
                                             parallel=False),
                              (1, 6))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=7,
                                             parallel=True),
                              (3, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3,
                                             max_procs=7,
                                             parallel=False),
                              (1, 7))


if __name__ == "__main__":
    ut.main()
