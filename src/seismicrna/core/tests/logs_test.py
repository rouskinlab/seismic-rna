import os
import unittest as ut
from itertools import product

from seismicrna.core.logs import (logger,
                                  Logger,
                                  Level,
                                  LoggerConfig,
                                  FILE_VERBOSITY,
                                  format_console_color,
                                  format_console_plain,
                                  format_logfile,
                                  exc_info,
                                  get_config,
                                  set_config,
                                  erase_config,
                                  restore_config)


class TestLoggerClass(ut.TestCase):

    def test_logger_class(self):
        self.assertIsInstance(logger, Logger)


class TestLevels(ut.TestCase):

    def test_levels(self):
        self.assertEqual(Level.FATAL, -3)
        self.assertEqual(Level.ERROR, -2)
        self.assertEqual(Level.WARNING, -1)
        self.assertEqual(Level.STATUS, 0)
        self.assertEqual(Level.PROCESS, 1)
        self.assertEqual(Level.ROUTINE, 2)
        self.assertEqual(Level.DETAIL, 3)


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


class TestExcInfo(ut.TestCase):

    @restore_config
    def test_exc_info(self):
        for verbosity in [Level.PROCESS,
                          Level.STATUS,
                          Level.WARNING,
                          Level.ERROR,
                          Level.FATAL]:
            set_config(verbosity=verbosity)
            self.assertFalse(exc_info())
        for verbosity in [Level.DETAIL,
                          Level.ROUTINE]:
            set_config(verbosity=verbosity)
            self.assertTrue(exc_info())


class TestLoggingRaiseOnError(ut.TestCase):

    @restore_config
    def test_raise_on_error(self):
        set_config(verbosity=(Level.FATAL - 1), raise_on_error=True)
        logger.warning("An error has occurred")
        self.assertRaisesRegex(RuntimeError,
                               "An error has occurred",
                               logger.error,
                               "An error has occurred")
        self.assertRaisesRegex(RuntimeError,
                               "A fatal error has occurred",
                               logger.fatal,
                               "A fatal error has occurred")
        logger.warning(ZeroDivisionError("Cannot divide by 0"))
        self.assertRaisesRegex(ZeroDivisionError,
                               "Cannot divide by 0",
                               logger.error,
                               ZeroDivisionError("Cannot divide by 0"))
        self.assertRaisesRegex(ZeroDivisionError,
                               "Cannot divide by 0",
                               logger.fatal,
                               ZeroDivisionError("Cannot divide by 0"))

    @restore_config
    def test_no_raise_on_error(self):
        set_config(verbosity=(Level.FATAL - 1), raise_on_error=False)
        # None of these calls should raise an error.
        logger.error("An error has occurred")
        logger.fatal("A fatal error has occurred")
        logger.error(ZeroDivisionError("Cannot divide by 0"))
        logger.fatal(ZeroDivisionError("Cannot divide by 0"))


class TestEraseConfig(ut.TestCase):

    def test_erase_config(self):
        set_config(verbosity=3, raise_on_error=True)
        self.assertEqual(logger.console_stream.filterer.verbosity, 3)
        self.assertTrue(logger.raise_on_error)
        erase_config()
        self.assertIsNone(logger.console_stream)
        self.assertIsNone(logger.file_stream)
        self.assertFalse(logger.raise_on_error)


class TestSetConfig(ut.TestCase):

    @restore_config
    def test_defaults(self):
        erase_config()
        self.assertIsNone(logger.console_stream)
        self.assertIsNone(logger.file_stream)
        self.assertFalse(logger.raise_on_error)
        set_config()
        self.assertFalse(logger.raise_on_error)
        self.assertEqual(logger.console_stream.filterer.verbosity, 0)
        self.assertIs(logger.console_stream.formatter.formatter,
                      format_console_color)

    @restore_config
    def test_verbosity(self):
        for verbosity in [-3, -2, -1, 0, 1, 2, 3]:
            set_config(verbosity=verbosity)
            self.assertEqual(logger.console_stream.filterer.verbosity,
                             verbosity)

    @restore_config
    def test_log_file(self):
        log_file_path = "test.log"
        msg1 = "Some logging text"
        msg2 = "More logging text"
        log_file = open(log_file_path, "x")
        try:
            set_config(log_file=log_file)
            self.assertEqual(logger.file_stream.filterer.verbosity,
                             FILE_VERBOSITY)
            self.assertIs(logger.file_stream.formatter.formatter,
                          format_logfile)
            self.assertIs(logger.file_stream.stream,
                          log_file)
            # Test logging a message to the file.
            logger.routine(msg1)
            logger.detail(msg2)
            log_file.close()
            with open(log_file_path) as log_file:
                lines = log_file.readlines()
                self.assertEqual(len(lines), 3)
                self.assertTrue(lines[0].startswith("LOGMSG> "))
                self.assertEqual(lines[1], f"{msg1}\n")
                self.assertEqual(lines[2], "\n")
        finally:
            log_file.close()
            try:
                os.remove(log_file_path)
            except FileNotFoundError:
                pass

    @restore_config
    def test_no_log_color(self):
        set_config(log_color=False)
        self.assertFalse(logger.raise_on_error)
        self.assertEqual(logger.console_stream.filterer.verbosity, 0)
        self.assertIs(logger.console_stream.formatter.formatter,
                      format_console_plain)

    @restore_config
    def test_raise_on_error(self):
        set_config(raise_on_error=True)
        self.assertTrue(logger.raise_on_error)


class TestGetConfig(ut.TestCase):

    @restore_config
    def test_get_config(self):
        log_file = os.path.abspath("test.log")
        options = {"verbosity": list(Level),
                   "log_file": [None, log_file],
                   "log_color": [True, False],
                   "raise_on_error": [False, True]}
        try:
            for choice in product(*options.values()):
                config = dict(zip(options.keys(), choice, strict=True))
                set_config(**config)
                result = get_config()
                self.assertIsInstance(result, LoggerConfig)
                self.assertDictEqual(result._asdict(), config)
        finally:
            try:
                os.remove(log_file)
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    ut.main()
