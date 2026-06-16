import unittest as ut
from unittest import mock

from seismicrna.core.logs import Level, logger, restore_config, set_config
from seismicrna.core.task import Task, calc_pool_size


class TestCalcPoolSize(ut.TestCase):
    def test_1_task(self):
        for num_cpus in range(1, 5):
            with self.subTest(num_cpus=num_cpus):
                expect = 1, num_cpus
                self.assertTupleEqual(calc_pool_size(1, num_cpus), expect)

    def test_1_proc(self):
        for n_tasks in range(1, 5):
            with self.subTest(n_tasks=n_tasks):
                expect = 1, 1
                self.assertTupleEqual(calc_pool_size(n_tasks, 1), expect)

    def test_multiple(self):
        self.assertTupleEqual(calc_pool_size(num_tasks=2, num_cpus=2), (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2, num_cpus=3), (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=2, num_cpus=4), (2, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=2, num_cpus=5), (2, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3, num_cpus=2), (2, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3, num_cpus=3), (3, 1))
        self.assertTupleEqual(calc_pool_size(num_tasks=3, num_cpus=6), (3, 2))
        self.assertTupleEqual(calc_pool_size(num_tasks=3, num_cpus=7), (3, 2))


class TestTaskDepth(ut.TestCase):
    @restore_config
    def test_child_starts_one_level_deeper(self):
        set_config(verbosity=Level.TRACE, log_color=False)
        observed = {}

        def record():
            observed["depth"] = len(logger.context_levels)
            return "done"

        # Create the task while the spawning process has 3 open contexts.
        logger.context_levels = [Level.INFO, Level.INFO, Level.INFO]
        task = Task(record)

        # Same process (PID matches): only begin() adds a context.
        self.assertEqual(task(), "done")
        same_process_depth = observed["depth"]

        # Child process (PID differs): the task reproduces the parent's
        # contexts plus one boundary level, then begin() adds one further.
        with mock.patch("seismicrna.core.task.getpid", return_value=task._pid + 1):
            self.assertEqual(task(), "done")
        child_depth = observed["depth"]

        self.assertEqual(child_depth, same_process_depth + 1)
        self.assertEqual(same_process_depth, 3)
        self.assertEqual(child_depth, 4)


if __name__ == "__main__":
    ut.main(verbosity=2)
