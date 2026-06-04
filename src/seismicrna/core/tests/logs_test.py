import io
import os
import unittest as ut
from itertools import product
from pathlib import Path
from unittest import mock

from seismicrna.core.logs import (
    logger,
    Logger,
    Level,
    LoggerConfig,
    Message,
    INDENT,
    FILE_VERBOSITY,
    format_console_color,
    format_console_plain,
    format_logfile,
    get_config,
    set_config,
    erase_config,
    restore_config,
)


class TestLoggerClass(ut.TestCase):
    def test_logger_class(self):
        self.assertIsInstance(logger, Logger)


class TestLevels(ut.TestCase):
    def test_levels(self):
        self.assertEqual(Level.ERROR, -2)
        self.assertEqual(Level.WARNING, -1)
        self.assertEqual(Level.INFO, 0)
        self.assertEqual(Level.DEBUG, 1)
        self.assertEqual(Level.TRACE, 2)
        self.assertEqual(
            [level.name for level in Level],
            ["ERROR", "WARNING", "INFO", "DEBUG", "TRACE"],
        )


class TestRestoreConfig(ut.TestCase):
    @restore_config
    def test_restore_config_no_error(self):
        def config_modifier_no_restore():
            set_config(verbosity=2)
            self.assertEqual(get_config().verbosity, 2)

        @restore_config
        def config_modifier_restore():
            set_config(verbosity=2)
            self.assertEqual(get_config().verbosity, 2)

        set_config(verbosity=-1)
        config_modifier_restore()
        self.assertEqual(get_config().verbosity, -1)
        config_modifier_no_restore()
        self.assertEqual(get_config().verbosity, 2)

    @restore_config
    def test_restore_config_error(self):
        def config_modifier_no_restore():
            set_config(verbosity=2)
            self.assertEqual(get_config().verbosity, 2)
            return 1 / 0

        @restore_config
        def config_modifier_restore():
            set_config(verbosity=2)
            self.assertEqual(get_config().verbosity, 2)
            return 1 / 0

        set_config(verbosity=-1)
        self.assertRaises(ZeroDivisionError, config_modifier_restore)
        self.assertEqual(get_config().verbosity, -1)
        self.assertRaises(ZeroDivisionError, config_modifier_no_restore)
        self.assertEqual(get_config().verbosity, 2)


class TestLoggingExitOnError(ut.TestCase):
    @restore_config
    def test_exit_on_error(self):
        set_config(verbosity=(Level.ERROR - 1), exit_on_error=True)
        # A message logged at the error level raises.
        self.assertRaisesRegex(
            RuntimeError, "An error has occurred", logger.error, "An error has occurred"
        )
        self.assertRaisesRegex(
            ZeroDivisionError,
            "Cannot divide by 0",
            logger.error,
            ZeroDivisionError("Cannot divide by 0"),
        )
        # A message logged below the error level never raises, even when its
        # content is an exception.
        logger.warning("An error has occurred")
        logger.warning(ZeroDivisionError("Cannot divide by 0"))

    @restore_config
    def test_no_exit_on_error(self):
        set_config(verbosity=(Level.ERROR - 1), exit_on_error=False)
        # None of these calls should raise an error.
        logger.error("An error has occurred")
        logger.error(ZeroDivisionError("Cannot divide by 0"))


class TestEraseConfig(ut.TestCase):
    def test_erase_config(self):
        set_config(verbosity=3, exit_on_error=False)
        self.assertEqual(logger.console_stream.filterer.verbosity, 3)
        self.assertFalse(logger.exit_on_error)
        erase_config()
        self.assertIsNone(logger.console_stream)
        self.assertIsNone(logger.file_stream)
        self.assertTrue(logger.exit_on_error)


class TestSetConfig(ut.TestCase):
    @restore_config
    def test_defaults(self):
        erase_config()
        self.assertIsNone(logger.console_stream)
        self.assertIsNone(logger.file_stream)
        self.assertTrue(logger.exit_on_error)
        set_config()
        self.assertTrue(logger.exit_on_error)
        self.assertEqual(logger.console_stream.filterer.verbosity, 0)
        self.assertIs(logger.console_stream.formatter.formatter, format_console_color)

    @restore_config
    def test_verbosity(self):
        for verbosity in [-3, -2, -1, 0, 1, 2, 3]:
            set_config(verbosity=verbosity)
            self.assertEqual(logger.console_stream.filterer.verbosity, verbosity)

    @restore_config
    def test_log_file(self):
        log_file_path = Path(os.path.abspath("test.log"))
        msg1 = "Some logging text"
        msg2 = "More logging text"
        try:
            set_config(log_file_path=log_file_path)
            self.assertEqual(logger.file_stream.filterer.verbosity, FILE_VERBOSITY)
            self.assertIs(logger.file_stream.formatter.formatter, format_logfile)
            self.assertEqual(logger.file_stream.file_path, log_file_path)
            # Test logging a message to the file.
            logger.debug(msg1)
            logger.trace(msg2)
            # The file stream must be closed explicitly to flush the log
            # messages to the file before reading with readlines().
            logger.file_stream.close()
            with open(log_file_path) as log_file:
                lines = log_file.readlines()
                self.assertEqual(len(lines), 6)
                self.assertTrue(
                    lines[0].startswith(f"LOGMSG> {Level.DEBUG.name} {os.getpid()}")
                )
                self.assertEqual(lines[1], f"{msg1}\n")
                self.assertEqual(lines[2], "\n")
                self.assertTrue(
                    lines[3].startswith(f"LOGMSG> {Level.TRACE.name} {os.getpid()}")
                )
                self.assertEqual(lines[4], f"{msg2}\n")
                self.assertEqual(lines[5], "\n")
        finally:
            try:
                os.remove(log_file_path)
            except FileNotFoundError:
                pass

    @restore_config
    def test_no_log_color(self):
        set_config(log_color=False)
        self.assertTrue(logger.exit_on_error)
        self.assertEqual(logger.console_stream.filterer.verbosity, 0)
        self.assertIs(logger.console_stream.formatter.formatter, format_console_plain)

    @restore_config
    def test_exit_on_error(self):
        set_config(exit_on_error=False)
        self.assertFalse(logger.exit_on_error)


class TestGetConfig(ut.TestCase):
    @restore_config
    def test_get_config(self):
        log_file_path = Path(os.path.abspath("test.log"))
        options = {
            "verbosity": list(Level),
            "log_file_path": [None, log_file_path],
            "log_color": [True, False],
            "exit_on_error": [False, True],
        }
        try:
            for choice in product(*options.values()):
                config = dict(zip(options.keys(), choice, strict=True))
                set_config(**config)
                result = get_config()
                self.assertIsInstance(result, LoggerConfig)
                self.assertDictEqual(result._asdict(), config)
        finally:
            try:
                os.remove(log_file_path)
            except FileNotFoundError:
                pass


class TestFormatConsole(ut.TestCase):
    def test_level_pid_message(self):
        message = Message(Level.INFO, "hello world")
        formatted = format_console_plain(message)
        self.assertIn(Level.INFO.name, formatted)
        self.assertIn(f"{os.getpid()}", formatted)
        self.assertTrue(formatted.endswith("hello world\n"))

    def test_depth_indents_message_only(self):
        # At depth 0 the message follows the prefix directly.
        flat = format_console_plain(Message(Level.TRACE, "msg"), 0)
        # At depth 2 the message is indented by two INDENT units, but the
        # level/PID prefix is unchanged.
        nested = format_console_plain(Message(Level.TRACE, "msg"), 2)
        prefix = f"{Level.TRACE.name: <8}{os.getpid()}  "
        self.assertEqual(flat, f"{prefix}msg\n")
        self.assertEqual(nested, f"{prefix}{INDENT * 2}msg\n")

    def test_multiline_message_aligned(self):
        # Continuation lines align under the message column (the level/PID
        # prefix width plus the depth indent), not under column 0.
        formatted = format_console_plain(Message(Level.TRACE, "line1\nline2"), 1)
        prefix_len = len(f"{Level.TRACE.name: <8}{os.getpid()}  ")
        out_lines = formatted.splitlines()
        self.assertTrue(out_lines[0].endswith("line1"))
        self.assertEqual(out_lines[1], " " * prefix_len + INDENT + "line2")

    def test_long_message_wrapped_to_terminal_width(self):
        # Force a known terminal width via the COLUMNS environment variable.
        with mock.patch.dict(os.environ, {"COLUMNS": "50"}):
            content = " ".join(["word"] * 40)
            formatted = format_console_plain(Message(Level.INFO, content))
        out_lines = formatted.splitlines()
        self.assertGreater(len(out_lines), 1)
        # Every line is wrapped to at most one less than the terminal width.
        for line in out_lines:
            self.assertLessEqual(len(line), 49)
        # No content is lost by wrapping.
        self.assertEqual(formatted.count("word"), 40)


class TestDepth(ut.TestCase):
    @restore_config
    def test_begin_pushes_and_restores(self):
        set_config(verbosity=Level.TRACE, log_color=False)
        self.assertEqual(len(logger.context_levels), 0)
        with logger.debug.begin("outer"):
            self.assertEqual(len(logger.context_levels), 1)
            with logger.debug.begin("inner"):
                self.assertEqual(len(logger.context_levels), 2)
            self.assertEqual(len(logger.context_levels), 1)
        self.assertEqual(len(logger.context_levels), 0)

    @restore_config
    def test_context_restored_on_exception(self):
        set_config(verbosity=Level.TRACE, log_color=False, exit_on_error=False)
        self.assertEqual(len(logger.context_levels), 0)
        with self.assertRaises(ValueError):
            with logger.debug.begin("task"):
                self.assertEqual(len(logger.context_levels), 1)
                raise ValueError("boom")
        self.assertEqual(len(logger.context_levels), 0)

    @restore_config
    def test_reconfigure_inside_begin_does_not_crash(self):
        # Reconfiguring logging while a context is open (as the `test` command
        # does) must not corrupt the stack or raise on exit.
        set_config(verbosity=Level.INFO, log_color=False)
        with logger.info.begin("outer"):
            self.assertEqual(len(logger.context_levels), 1)
            set_config(verbosity=Level.DEBUG, log_color=False)
            self.assertEqual(len(logger.context_levels), 1)
        self.assertEqual(len(logger.context_levels), 0)

    def test_depth_counts_only_visible_contexts(self):
        # A message enclosed by an INFO context and a DEBUG context.
        message = Message(Level.INFO, "x", context_levels=[Level.INFO, Level.DEBUG])
        # At INFO verbosity the DEBUG context is hidden, so it adds no indent.
        self.assertEqual(message.depth(Level.INFO), 1)
        # At DEBUG verbosity both contexts are visible and both indent.
        self.assertEqual(message.depth(Level.DEBUG), 2)
        # At WARNING verbosity neither context is visible.
        self.assertEqual(message.depth(Level.WARNING), 0)

    @restore_config
    def test_hidden_context_does_not_indent_console(self):
        # At INFO verbosity, an INFO message inside a DEBUG context is not
        # indented, because the enclosing DEBUG "Began" line is itself hidden.
        set_config(verbosity=Level.INFO, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            with logger.debug.begin("hidden outer"):
                logger.info("visible message")
        out = stream.getvalue()
        self.assertNotIn("hidden outer", out)
        prefix = f"{Level.INFO.name: <8}{os.getpid()}  "
        self.assertIn(f"{prefix}visible message", out)


class TestBeginEnd(ut.TestCase):
    @restore_config
    def test_logs_began_and_ended(self):
        set_config(verbosity=Level.TRACE, log_color=False)
        stream = io.StringIO()
        # logs.py binds `stderr` at import, so patch it directly to capture
        # console output (redirect_stderr would not affect it).
        with mock.patch("seismicrna.core.logs.stderr", stream):
            with logger.debug.begin("doing x"):
                logger.trace("progress")
        output = stream.getvalue()
        self.assertIn("Began doing x", output)
        self.assertIn("Ended doing x", output)
        # The nested message is indented by one level.
        self.assertIn(f"{INDENT}progress", output)

    @restore_config
    def test_ended_skipped_on_exception(self):
        set_config(verbosity=Level.TRACE, log_color=False, exit_on_error=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            with self.assertRaises(ValueError):
                with logger.debug.begin("doing x"):
                    raise ValueError("boom")
        output = stream.getvalue()
        self.assertIn("Began doing x", output)
        self.assertNotIn("Ended doing x", output)


class TestVerbatimCall(ut.TestCase):
    @restore_config
    def test_verbatim_call_still_works(self):
        set_config(verbosity=Level.TRACE, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            logger.trace("just a message")
        self.assertIn("just a message", stream.getvalue())
        self.assertIsInstance(logger, Logger)


if __name__ == "__main__":
    ut.main(verbosity=2)
