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
    INDENT,
    format_body,
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
        self.assertRaisesRegex(
            RuntimeError, "An error has occurred", logger.error, "An error has occurred"
        )
        self.assertRaisesRegex(
            ZeroDivisionError,
            "Cannot divide by 0",
            logger.error,
            ZeroDivisionError("Cannot divide by 0"),
        )
        # Messages below the error level never raise, even with exception content.
        logger.warning("An error has occurred")
        logger.warning(ZeroDivisionError("Cannot divide by 0"))

    @restore_config
    def test_no_exit_on_error(self):
        set_config(verbosity=(Level.ERROR - 1), exit_on_error=False)
        logger.error("An error has occurred")
        logger.error(ZeroDivisionError("Cannot divide by 0"))


class TestEraseConfig(ut.TestCase):
    def test_erase_config(self):
        set_config(verbosity=3, exit_on_error=False)
        self.assertEqual(logger.verbosity, 3)
        self.assertFalse(logger.exit_on_error)
        erase_config()
        self.assertFalse(logger.console_enabled)
        self.assertIsNone(logger.file_stream)
        self.assertTrue(logger.exit_on_error)


class TestSetConfig(ut.TestCase):
    @restore_config
    def test_defaults(self):
        erase_config()
        self.assertFalse(logger.console_enabled)
        self.assertIsNone(logger.file_stream)
        self.assertTrue(logger.exit_on_error)
        set_config()
        self.assertTrue(logger.exit_on_error)
        self.assertEqual(logger.verbosity, 0)
        self.assertTrue(logger.log_color)

    @restore_config
    def test_verbosity(self):
        for verbosity in [-3, -2, -1, 0, 1, 2, 3]:
            set_config(verbosity=verbosity)
            self.assertEqual(logger.verbosity, verbosity)

    @restore_config
    def test_log_file_same_verbosity(self):
        """File records at the same verbosity as the console."""
        log_file_path = Path(os.path.abspath("test_same_verbosity.log"))
        msg_info = "An info message"
        msg_debug = "Some debug text"
        msg_trace = "Some trace text"
        try:
            # Set verbosity to INFO — DEBUG and TRACE should NOT appear in file.
            set_config(
                verbosity=Level.INFO, log_file_path=log_file_path, log_color=False
            )
            stream = io.StringIO()
            with mock.patch("seismicrna.core.logs.stderr", stream):
                logger.info(msg_info)
                logger.debug(msg_debug)
                logger.trace(msg_trace)
            logger.file_stream.close()
            with open(log_file_path) as log_file:
                content = log_file.read()
            self.assertIn(msg_info, content)
            self.assertNotIn(msg_debug, content)
            self.assertNotIn(msg_trace, content)
        finally:
            try:
                os.remove(log_file_path)
            except FileNotFoundError:
                pass

    @restore_config
    def test_log_file_timestamp_prefix(self):
        """File lines carry a timestamp prefix; body matches console body."""
        log_file_path = Path(os.path.abspath("test.log"))
        msg = "Hello log file"
        # Timestamp format is always "YYYY-%m-%d %H:%M:%S" (19 chars) + "  " = 28.
        ts_prefix_len = len("YYYY-MM-DD HH:MM:SS.ffffff  ")
        try:
            set_config(
                verbosity=Level.DEBUG, log_file_path=log_file_path, log_color=False
            )
            console_stream = io.StringIO()
            with mock.patch("seismicrna.core.logs.stderr", console_stream):
                logger.debug(msg)
            logger.file_stream.close()
            with open(log_file_path) as log_file:
                file_content = log_file.read()
            console_body = console_stream.getvalue()
            # File line = "<YYYY-MM-DD HH:MM:SS.ffffff> " + console body
            self.assertTrue(file_content.startswith("20"))  # starts with year
            # The body portion (after the timestamp prefix) matches the console.
            self.assertEqual(file_content[ts_prefix_len:], console_body)
        finally:
            try:
                os.remove(log_file_path)
            except FileNotFoundError:
                pass

    @restore_config
    def test_no_log_color(self):
        set_config(log_color=False)
        self.assertFalse(logger.log_color)

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


class TestFormatBody(ut.TestCase):
    def test_level_pid_message(self):
        body = format_body(Level.INFO, "hello world", 0)
        self.assertIn(Level.INFO.name.capitalize(), body)
        self.assertIn(f"{os.getpid()}", body)
        self.assertTrue(body.endswith("hello world\n"))

    def test_depth_indents_message_only(self):
        flat = format_body(Level.TRACE, "msg", 0)
        nested = format_body(Level.TRACE, "msg", 2)
        prefix = f"{Level.TRACE.name.capitalize():<7} {str(os.getpid())[:7]:>7} "
        self.assertEqual(flat, f"{prefix}msg\n")
        self.assertEqual(nested, f"{prefix}{INDENT * 2}msg\n")

    def test_multiline_message_aligned(self):
        formatted = format_body(Level.TRACE, "line1\nline2", 1)
        prefix_len = len(
            f"{Level.TRACE.name.capitalize():<7} {str(os.getpid())[:7]:>7} "
        )
        out_lines = formatted.splitlines()
        self.assertTrue(out_lines[0].endswith("line1"))
        self.assertEqual(out_lines[1], " " * prefix_len + INDENT + "line2")

    def test_long_message_wrapped_to_terminal_width(self):
        with mock.patch.dict(os.environ, {"COLUMNS": "50"}):
            body = format_body(Level.INFO, " ".join(["word"] * 40), 0)
        out_lines = body.splitlines()
        self.assertGreater(len(out_lines), 1)
        for line in out_lines:
            self.assertLessEqual(len(line), 49)
        self.assertEqual(body.count("word"), 40)


class TestLazyFormatting(ut.TestCase):
    """Template is formatted only when the level is visible."""

    @restore_config
    def test_skips_format_when_level_hidden(self):
        set_config(verbosity=Level.INFO, log_color=False, exit_on_error=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            # This would raise KeyError if .format() were called — field "bad"
            # is not in the positional args.  At TRACE verbosity (hidden at INFO)
            # it must be silently skipped.
            logger.trace("value = {bad}", "oops")
        self.assertEqual(stream.getvalue(), "")

    @restore_config
    def test_formats_when_level_visible(self):
        set_config(verbosity=Level.TRACE, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            logger.trace("x = {}, y = {}", 1, 2)
        self.assertIn("x = 1, y = 2", stream.getvalue())

    @restore_config
    def test_keyword_args(self):
        set_config(verbosity=Level.DEBUG, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            logger.debug("ref={ref} reg={reg}", ref="myref", reg="myreg")
        out = stream.getvalue()
        self.assertIn("ref=myref", out)
        self.assertIn("reg=myreg", out)

    @restore_config
    def test_plain_string_no_args(self):
        set_config(verbosity=Level.INFO, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            logger.info("plain message with no args")
        self.assertIn("plain message with no args", stream.getvalue())

    @restore_config
    def test_exception_object(self):
        set_config(verbosity=Level.ERROR, log_color=False, exit_on_error=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            logger.error(ValueError("something broke"))
        self.assertIn("something broke", stream.getvalue())


class TestDepth(ut.TestCase):
    @restore_config
    def test_begin_pushes_and_restores(self):
        set_config(verbosity=Level.TRACE, log_color=False)
        initial = len(logger.context_levels)
        with logger.debug.single_context("outer"):
            self.assertEqual(len(logger.context_levels), initial + 1)
            with logger.debug.single_context("inner"):
                self.assertEqual(len(logger.context_levels), initial + 2)
            self.assertEqual(len(logger.context_levels), initial + 1)
        self.assertEqual(len(logger.context_levels), initial)

    @restore_config
    def test_context_restored_on_exception(self):
        set_config(verbosity=Level.TRACE, log_color=False, exit_on_error=False)
        initial = len(logger.context_levels)
        with self.assertRaises(ValueError):
            with logger.debug.single_context("task"):
                self.assertEqual(len(logger.context_levels), initial + 1)
                raise ValueError("boom")
        self.assertEqual(len(logger.context_levels), initial)

    @restore_config
    def test_reconfigure_inside_begin_does_not_crash(self):
        set_config(verbosity=Level.INFO, log_color=False)
        initial = len(logger.context_levels)
        with logger.info.single_context("outer"):
            self.assertEqual(len(logger.context_levels), initial + 1)
            set_config(verbosity=Level.DEBUG, log_color=False)
            self.assertEqual(len(logger.context_levels), initial + 1)
        self.assertEqual(len(logger.context_levels), initial)

    @restore_config
    def test_hidden_context_does_not_indent_console(self):
        # At INFO verbosity, an INFO message inside a DEBUG context is not
        # indented more than one outside it, because the enclosing DEBUG
        # context is itself hidden.
        set_config(verbosity=Level.INFO, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            logger.info("baseline")
            with logger.debug.single_context("hidden outer"):
                logger.info("visible message")
        out = stream.getvalue()
        self.assertNotIn("hidden outer", out)
        lines = out.splitlines()
        # Both lines must end the same way (same indentation from outer context).
        self.assertEqual(
            lines[0].split("baseline")[0], lines[1].split("visible message")[0]
        )

    @restore_config
    def test_visible_context_indents_nested_message(self):
        # At DEBUG verbosity, messages inside a DEBUG context are indented.
        set_config(verbosity=Level.DEBUG, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            with logger.debug.single_context("outer"):
                logger.debug("inner message")
        out = stream.getvalue()
        # The "inner message" line must contain one level of indentation.
        inner_line = [line for line in out.splitlines() if "inner message" in line][0]
        self.assertIn(INDENT, inner_line)


class TestBeginEnd(ut.TestCase):
    @restore_config
    def test_logs_began_and_ended(self):
        set_config(verbosity=Level.TRACE, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            with logger.debug.double_context("doing x"):
                logger.trace("progress")
        output = stream.getvalue()
        self.assertIn("Began doing x", output)
        self.assertIn("Ended doing x", output)
        self.assertIn(f"{INDENT}progress", output)

    @restore_config
    def test_ended_skipped_on_exception(self):
        set_config(verbosity=Level.TRACE, log_color=False, exit_on_error=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            with self.assertRaises(ValueError):
                with logger.debug.double_context("doing x"):
                    raise ValueError("boom")
        output = stream.getvalue()
        self.assertIn("Began doing x", output)
        self.assertNotIn("Ended doing x", output)

    @restore_config
    def test_begin_with_format_args(self):
        set_config(verbosity=Level.DEBUG, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            with logger.debug.double_context("processing {} items", 42):
                pass
        out = stream.getvalue()
        self.assertIn("Began processing 42 items", out)
        self.assertIn("Ended processing 42 items", out)

    @restore_config
    def test_begin_hidden_does_not_format(self):
        # begin() at DEBUG level is hidden at INFO verbosity — template must
        # not be formatted.
        set_config(verbosity=Level.INFO, log_color=False)
        stream = io.StringIO()
        with mock.patch("seismicrna.core.logs.stderr", stream):
            # Passing bad positional arg — would raise if .format() were called.
            with logger.debug.double_context("value = {bad}", "oops"):
                logger.info("still visible")
        out = stream.getvalue()
        self.assertNotIn("value", out)
        self.assertIn("still visible", out)


if __name__ == "__main__":
    ut.main(verbosity=2)
