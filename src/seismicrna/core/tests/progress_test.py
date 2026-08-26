import io
import unittest as ut
from unittest import mock

from seismicrna.core import logs, progress
from seismicrna.core.logs import Level
from seismicrna.core.logs import restore_config as restore_log_config
from seismicrna.core.logs import set_config as set_log_config
from seismicrna.core.progress import (
    ConsoleStream,
    DEFAULT_ENABLED,
    DESC_CUT,
    DESC_WIDTH,
    ERASE_BARS,
    ProgressBar,
    ProgressConfig,
    close_bars,
    erase_config,
    format_desc,
    get_config,
    restore_config,
    set_config,
    set_shared_config,
    write_console,
)


class Terminal(io.StringIO):
    """Stream that claims to be a terminal, as stderr is in a console."""

    def isatty(self):
        return True


class ProgressStateTestCase(ut.TestCase):
    """Test case that leaves the progress bars as it found them.

    Progress bars are global to the process, so a test that reconfigures them
    would otherwise close and disable the progress bar of the test command
    that is running it. Each test starts with no progress bars configured and
    ends with the previous configuration (and any bar that it was drawing)
    back in place."""

    def setUp(self):
        self._state = (
            progress._enabled,
            progress._pid,
            progress._shared,
            list(progress._open_bars),
            logs.console_hook,
        )
        progress._enabled = False
        progress._shared = False
        progress._open_bars.clear()
        logs.console_hook = None

    def tearDown(self):
        # Close only the bars that this test itself opened, and take them
        # off the console rather than leaving them for the next test.
        for bar in progress._live_bars():
            bar.close()
        for bar in list(progress._open_bars):
            bar._erase()
        progress._open_bars.clear()
        (
            progress._enabled,
            progress._pid,
            progress._shared,
            open_bars,
            logs.console_hook,
        ) = self._state
        progress._open_bars[:] = open_bars

    def enable_bars(self):
        """Draw progress bars on a fake terminal, and return that terminal."""
        stream = Terminal()
        patcher = mock.patch.object(logs, "stderr", stream)
        patcher.start()
        # Stop patching only after tearDown has closed this test's bars, so
        # that closing them writes to the fake terminal, not the real one.
        self.addCleanup(patcher.stop)
        set_config(True, Level.INFO)
        return stream


class TestFormatDesc(ut.TestCase):
    def test_padded_to_a_fixed_width(self):
        # Every description is the same width, so that the gauges after them
        # line up in a column and never move as tasks begin and end.
        for label in ["", "aligning", "identifying mutations"]:
            with self.subTest(label=label):
                self.assertEqual(len(format_desc(label)), DESC_WIDTH)
                self.assertTrue(format_desc(label).startswith(label))

    def test_long_label_cut_short(self):
        label = "x" * (DESC_WIDTH + 10)
        desc = format_desc(label)
        self.assertEqual(len(desc), DESC_WIDTH)
        self.assertTrue(desc.endswith(DESC_CUT))


class TestSetConfig(ProgressStateTestCase):
    def test_default_disabled(self):
        self.assertFalse(DEFAULT_ENABLED)
        self.assertEqual(get_config(), ProgressConfig(enabled=False))
        self.assertIsNone(logs.console_hook)

    def test_enable_on_terminal(self):
        with mock.patch.object(logs, "stderr", Terminal()):
            set_config(True, Level.INFO)
            self.assertTrue(get_config().enabled)
            self.assertIs(logs.console_hook, write_console)

    def test_disable_off_terminal(self):
        # io.StringIO is not a terminal, so bars would be meaningless.
        with mock.patch.object(logs, "stderr", io.StringIO()):
            set_config(True, Level.INFO)
            self.assertFalse(get_config().enabled)
            self.assertIsNone(logs.console_hook)

    def test_disable_low_verbosity(self):
        with mock.patch.object(logs, "stderr", Terminal()):
            for verbosity in [Level.WARNING, Level.ERROR]:
                with self.subTest(verbosity=verbosity):
                    set_config(True, verbosity)
                    self.assertFalse(get_config().enabled)
            for verbosity in [Level.INFO, Level.DEBUG, Level.TRACE]:
                with self.subTest(verbosity=verbosity):
                    set_config(True, verbosity)
                    self.assertTrue(get_config().enabled)

    def test_erase_config_removes_hook(self):
        with mock.patch.object(logs, "stderr", Terminal()):
            set_config(True, Level.INFO)
            erase_config()
            self.assertFalse(get_config().enabled)
            self.assertIsNone(logs.console_hook)

    def test_erase_config_leaves_the_bars_on_the_console(self):
        stream = Terminal()
        with mock.patch.object(logs, "stderr", stream):
            set_config(True, Level.INFO)
            bar = ProgressBar("task", 1)
            bar.tick()
            stream.truncate(0)
            stream.seek(0)
            erase_config()
        # The bars stay on the console as a record of the run; the cursor is
        # put below them so that nothing written next overwrites them.
        self.assertEqual(progress._open_bars, [])
        output = stream.getvalue()
        self.assertIn(format_desc("task"), output)
        self.assertTrue(output.endswith("\n"))

    def test_restore_config(self):
        with mock.patch.object(logs, "stderr", Terminal()):

            @restore_config
            def disable_bars():
                erase_config()
                self.assertFalse(get_config().enabled)

            set_config(True, Level.INFO)
            disable_bars()
            self.assertTrue(get_config().enabled)
            self.assertIs(logs.console_hook, write_console)


class TestWriteConsole(ProgressStateTestCase):
    def test_write_without_bars(self):
        stream = io.StringIO()
        with mock.patch.object(logs, "stderr", stream):
            write_console("message\n")
        self.assertEqual(stream.getvalue(), "message\n")

    def test_write_from_other_process(self):
        # A process that shares the console must erase the bars that another
        # process drew before writing over them.
        stream = Terminal()
        with mock.patch.object(logs, "stderr", stream):
            set_shared_config(True)
            self.assertIs(logs.console_hook, write_console)
            write_console("message\n")
        self.assertEqual(stream.getvalue(), f"{ERASE_BARS}message\n")

    def test_child_process_does_not_draw(self):
        stream = Terminal()
        with mock.patch.object(logs, "stderr", stream):
            set_config(True, Level.INFO)
            self.assertTrue(progress._owns_console())
            # Pretend that this process was forked from the one that
            # enabled the progress bars.
            with mock.patch.object(progress, "getpid", lambda: progress._pid + 1):
                self.assertFalse(progress._owns_console())
                stream.truncate(0)
                stream.seek(0)
                write_console("message\n")
                # The bar belongs to the process that this one was forked
                # from, so this process erases it rather than drawing it,
                # even before its task has told it that it shares the
                # console.
                self.assertTrue(progress._shares_console())
                self.assertEqual(stream.getvalue(), f"{ERASE_BARS}message\n")

    def test_failure_disables_bars(self):
        stream = Terminal()
        with mock.patch.object(logs, "stderr", stream):
            set_config(True, Level.INFO)
            with mock.patch.object(
                progress.tqdm, "write", side_effect=OSError("broken")
            ):
                write_console("message\n")
            # Logging must continue even though drawing failed.
            self.assertEqual(stream.getvalue(), "message\n")
            self.assertFalse(get_config().enabled)
            self.assertIsNone(logs.console_hook)


class TestProgressBar(ProgressStateTestCase):
    def test_disabled_bar_does_nothing(self):
        stream = io.StringIO()
        with mock.patch.object(logs, "stderr", stream):
            with ProgressBar("task", 3) as bar:
                bar.tick()
                bar.tick(2)
                self.assertEqual(bar.count, 3)
        self.assertEqual(stream.getvalue(), "")

    def test_counts_are_clamped(self):
        with ProgressBar("task", 2) as bar:
            bar.tick(5)
            self.assertEqual(bar.count, 2)
            bar.update(1)
            self.assertEqual(bar.count, 1)

    def test_set_label(self):
        with ProgressBar("before", 1) as bar:
            bar.set_label("after")
            self.assertEqual(bar.label, "after")
            bar.tick(label="last")
            self.assertEqual(bar.label, "last")

    def test_numbers_count_from_one(self):
        # The number shown is the item being worked on, not the number of
        # items finished, so it starts at 1 rather than 0.
        self.enable_bars()
        with ProgressBar("task", 3) as bar:
            self.assertEqual(bar.number, 1)
            bar.tick()
            self.assertEqual(bar.number, 2)
            bar.tick()
            self.assertEqual(bar.number, 3)
            # After the last item finishes, no further item begins.
            bar.tick()
            self.assertEqual(bar.count, 3)
            self.assertEqual(bar.number, 3)

    def test_no_items_numbered_zero(self):
        with ProgressBar("task", 0) as bar:
            self.assertEqual(bar.number, 0)

    def test_close_is_idempotent(self):
        self.enable_bars()
        bar = ProgressBar("task", 1)
        bar.close()
        bar.close()
        self.assertTrue(bar.done)
        self.assertEqual(progress._open_bars, [bar])

    def test_nested_bars_close_out_of_order(self):
        # A bar opened inside a generator can outlive the bar that was
        # opened after it, so bars do not always close in reverse order.
        self.enable_bars()
        outer = ProgressBar("outer", 1)
        inner = ProgressBar("inner", 1)
        self.assertEqual(progress._open_bars, [outer, inner])
        outer.close()
        self.assertTrue(outer.done)
        self.assertFalse(inner.done)
        # Both keep their lines, in the order their tasks began.
        self.assertEqual(progress._open_bars, [outer, inner])
        inner.close()
        self.assertTrue(inner.done)
        self.assertEqual(progress._open_bars, [outer, inner])

    def test_nested_tasks_get_their_own_lines(self):
        # A step of the workflow stays visible on its own line above the
        # sample it is working on.
        self.enable_bars()
        with ProgressBar("wf", 5) as step:
            step.tick(label="align")
            with ProgressBar("aligning reads", 2) as inner:
                self.assertTrue(step.shown)
                self.assertTrue(inner.shown)
                # Each bar names only itself, since each has its own line.
                self.assertEqual(step._format_desc(), format_desc("align"))
                self.assertEqual(inner._format_desc(), format_desc("aligning reads"))
                # The step is above the sample it is working on.
                self.assertLess(step._bar.pos, inner._bar.pos)

    def test_full_block_sends_finished_bars_to_the_log(self):
        # The block may not fill the console, or there would be no room for
        # the log, so once it is full each further task goes to the log.
        stream = self.enable_bars()
        with mock.patch.object(progress, "_max_bars", lambda: 1):
            with ProgressBar("wf", 5):
                stream.truncate(0)
                stream.seek(0)
                with ProgressBar("aligning reads", 2) as bar:
                    bar.tick(2)
                # The block was already full, so this task gave up its line.
                self.assertNotIn(bar, progress._open_bars)
                self.assertIn(format_desc("aligning reads"), stream.getvalue())

    def test_warning_written_above_the_bar(self):
        # Python writes warnings straight to standard error, bypassing the
        # logger, so they would otherwise land on the bar.
        import warnings

        stream = self.enable_bars()
        with ProgressBar("task", 2):
            stream.truncate(0)
            stream.seek(0)
            warnings.warn("something odd", RuntimeWarning)
            output = stream.getvalue()
        i = output.index("something odd")
        # Nothing of the bar precedes the warning on its line: the bar was
        # cleared first, rather than the warning landing on the end of it.
        line_start = output[:i].rsplit("\r", 1)[-1].rsplit("\n", 1)[-1]
        self.assertNotIn("|", line_start)
        # The bar is drawn again below the warning.
        self.assertIn(format_desc("task"), output[i:])

    def test_refresh_draws_every_bar(self):
        # A task in another process erases the whole block of bars before
        # writing a message, so drawing only the innermost bar again would
        # leave the tasks that contain it blank.
        stream = self.enable_bars()
        with ProgressBar("outer", 2), ProgressBar("inner", 2):
            stream.truncate(0)
            stream.seek(0)
            progress.refresh_all()
            output = stream.getvalue()
            self.assertIn(format_desc("outer"), output)
            self.assertIn(format_desc("inner"), output)

    def test_finished_bars_stay_and_the_stack_grows_downward(self):
        # Every task that finishes keeps its line, and the next task takes
        # the line below it, so the block grows downward in the order the
        # tasks began and ends up showing the whole run.
        self.enable_bars()
        with ProgressBar("wf", 5) as step:
            with ProgressBar("aligning reads", 1) as first:
                first.tick()
            self.assertTrue(first.done)
            with ProgressBar("identifying mutations", 2) as second:
                # The task that contains them is still on top, the finished
                # task keeps its line, and the new task is below both.
                self.assertEqual(progress._open_bars, [step, first, second])
                self.assertLess(step._bar.pos, first._bar.pos)
                self.assertLess(first._bar.pos, second._bar.pos)

    def test_finished_bars_survive_the_end_of_a_command(self):
        # The record of a task that finished must outlive the command that
        # ran it, so that the console shows the whole run at the end.
        self.enable_bars()

        @close_bars
        def command():
            with ProgressBar("aligning reads", 1) as bar:
                bar.tick()

        command()
        self.assertEqual(len(progress._open_bars), 1)
        self.assertTrue(progress._open_bars[0].done)

    def test_unnamed_task_gets_no_bar(self):
        # The parts of one sample's analysis are not named, so they neither
        # draw a bar nor take the line from the task that contains them.
        stream = self.enable_bars()
        with ProgressBar("named", 2) as named:
            stream.truncate(0)
            stream.seek(0)
            with ProgressBar("", 5) as unnamed:
                unnamed.tick()
                self.assertEqual(unnamed.count, 1)
                self.assertFalse(unnamed.shown)
                # The named task still has the only line.
                self.assertEqual(progress._open_bars, [named])
            self.assertEqual(stream.getvalue(), "")


class TestCloseBars(ProgressStateTestCase):
    def test_closes_bars_left_open(self):
        self.enable_bars()

        @close_bars
        def leak():
            ProgressBar("leaked", 1)
            return "done"

        outer = ProgressBar("kept", 1)
        self.assertEqual(leak(), "done")
        # The task the function left running is finished, but its bar stays
        # on the console; the task that contains it is still running.
        self.assertEqual(len(progress._open_bars), 2)
        self.assertFalse(outer.done)
        self.assertTrue(progress._open_bars[1].done)

    def test_closes_bars_after_error(self):
        self.enable_bars()

        @close_bars
        def fail():
            ProgressBar("leaked", 1)
            raise ValueError("failed")

        self.assertRaises(ValueError, fail)
        self.assertEqual(len(progress._open_bars), 1)
        self.assertTrue(progress._open_bars[0].done)


class TestConsoleStream(ProgressStateTestCase):
    def test_passthrough_without_bars(self):
        stream = io.StringIO()
        with mock.patch.object(logs, "stderr", stream):
            with ConsoleStream() as out:
                out.write("par")
                out.flush()
                # With no bars to protect, partial writes pass straight
                # through, so that output appears as soon as it is written.
                self.assertEqual(stream.getvalue(), "par")
                out.write("tial\n")
        self.assertEqual(stream.getvalue(), "partial\n")

    def test_buffers_partial_lines(self):
        stream = Terminal()
        with mock.patch.object(logs, "stderr", stream):
            set_config(True, Level.INFO)
            with ConsoleStream() as out:
                out.write("..")
                out.flush()
                # The line is incomplete, so it is held back: writing it now
                # would leave the cursor where the bars are drawn.
                self.assertEqual(stream.getvalue(), "")
                out.write("..\nOK\n")
                self.assertEqual(stream.getvalue(), "....\nOK\n")
                out.write("trailing")
            # Closing the stream writes whatever is left.
            self.assertEqual(stream.getvalue(), "....\nOK\ntrailing\n")


class TestLoggingWithBars(ProgressStateTestCase):
    @restore_log_config
    def test_log_message_goes_above_the_bar(self):
        stream = Terminal()
        with mock.patch.object(logs, "stderr", stream):
            set_log_config(Level.INFO, None, False, False)
            set_config(True, Level.INFO)
            with ProgressBar("task", 1) as bar:
                logs.logger.info("hello")
                bar.tick()
        output = stream.getvalue()
        # The bar was drawn, and the message was written at the start of a
        # line (after a carriage return that erased the bar) rather than
        # onto the end of the line where the bar is drawn.
        self.assertIn(format_desc("task"), output)
        self.assertIn("\rInfo", output)
        self.assertIn("hello\n", output)
        # The bar stays on the console when the configuration is erased,
        # with the cursor moved below it.
        with mock.patch.object(logs, "stderr", stream):
            stream.truncate(0)
            stream.seek(0)
            erase_config()
        self.assertTrue(stream.getvalue().endswith("\n"))


if __name__ == "__main__":
    ut.main(verbosity=2)
