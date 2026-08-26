import unittest as ut
from concurrent.futures import Future
from unittest import mock

from seismicrna.core import logs, progress
from seismicrna.core import task as task_module
from seismicrna.core.logs import Level, logger, restore_config, set_config
from seismicrna.core.progress import ProgressBar
from seismicrna.core.progress import erase_config as erase_progress_config
from seismicrna.core.progress import set_config as set_progress_config
from seismicrna.core.tests.progress_test import ProgressStateTestCase, Terminal
from seismicrna.core.task import (
    Task,
    _iter_finished,
    as_list_of_tuples,
    calc_pool_size,
    dispatch,
)


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


class TestTaskDepth(ProgressStateTestCase):
    """Test the logging context of a task in a child process.

    This derives from ProgressStateTestCase because pretending to be a child
    process also configures this process to share its console, which would
    otherwise stop the test command from drawing its own progress bar."""

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


class TestIterFinished(ProgressStateTestCase):
    """Test the loop that yields finished tasks and redraws the bar."""

    @staticmethod
    def _futures(values: list):
        futures = list()
        for value in values:
            future = Future()
            future.set_result(value)
            futures.append(future)
        return futures

    def test_ordered(self):
        futures = self._futures([1, 2, 3])
        bar = ProgressBar("test", 3)
        results = [f.result() for f in _iter_finished(futures, True, bar)]
        self.assertListEqual(results, [1, 2, 3])
        bar.close()

    def test_unordered(self):
        futures = self._futures([1, 2, 3])
        bar = ProgressBar("test", 3)
        results = [f.result() for f in _iter_finished(futures, False, bar)]
        self.assertCountEqual(results, [1, 2, 3])
        bar.close()

    def test_unnamed_task_still_draws_the_bars(self):
        # A task with no bar of its own still runs while the tasks that
        # contain it have bars on the console, and its own tasks write
        # messages over them, so those bars must still be drawn again.
        self.enable_bars()
        with ProgressBar("outer", 2):
            unnamed = ProgressBar("", 2)
            self.assertFalse(unnamed.shown)
            drawn = list()
            with mock.patch.object(task_module, "refresh_all", lambda: drawn.append(1)):
                got = [
                    f.result()
                    for f in _iter_finished(self._futures([1]), False, unnamed)
                ]
            self.assertListEqual(got, [1])
            self.assertTrue(drawn)

    def test_refreshes_while_waiting(self):
        # A future that never finishes must not stop the bar from being
        # redrawn, since a task in another process may have erased it.
        futures = self._futures([1]) + [Future()]
        bar = mock.MagicMock()
        finished = list()
        for future in _iter_finished(futures, False, bar):
            finished.append(future.result())
            if len(finished) == 1:
                break
        self.assertListEqual(finished, [1])


def _square(x: int):
    """Square a number (defined at module level so it can be pickled)."""
    return x * x


class TestDispatchProgress(ProgressStateTestCase):
    """Test that progress bars do not change the results of dispatch."""

    def _dispatch_squares(self, num_cpus: int):
        return dispatch(
            _square,
            num_cpus=num_cpus,
            pass_num_cpus=False,
            as_list=True,
            ordered=True,
            raise_on_error=True,
            label="squaring",
            args=as_list_of_tuples(range(6)),
        )

    def test_same_results_with_and_without_bars(self):
        expect = [x * x for x in range(6)]
        for num_cpus in [1, 4]:
            with self.subTest(num_cpus=num_cpus):
                erase_progress_config()
                self.assertListEqual(self._dispatch_squares(num_cpus), expect)
                with mock.patch.object(logs, "stderr", Terminal()):
                    set_progress_config(True, Level.INFO)
                    self.assertListEqual(self._dispatch_squares(num_cpus), expect)
                # Leave no bar drawing on the real console for the next case.
                erase_progress_config()

    def test_bar_is_finished_afterwards(self):
        # The bar keeps its line as a record of the task, but the task is
        # marked finished so that nothing tries to advance it again.
        with mock.patch.object(logs, "stderr", Terminal()):
            set_progress_config(True, Level.INFO)
            self._dispatch_squares(4)
            self.assertEqual(len(progress._open_bars), 1)
            self.assertTrue(progress._open_bars[0].done)
            self.assertEqual(progress._live_bars(), [])


if __name__ == "__main__":
    ut.main(verbosity=2)
