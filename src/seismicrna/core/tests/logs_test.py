import logging
import os
import unittest as ut
from itertools import product

from seismicrna.core.logs import (ColorFormatter,
                                  LoggerConfig,
                                  RaisableLogger,
                                  exc_info,
                                  get_top_logger,
                                  get_verbosity,
                                  get_config,
                                  set_config,
                                  erase_config,
                                  restore_config)
from seismicrna.core.tests import __name__ as __parent_name__

# Use top_logger any time an ATTRIBUTE is being checked, e.g. the list
# of handlers or the verbosity level.
top_logger = get_top_logger()
# Use test_logger any time a BEHAVIOR is being checked, e.g. the result
# of logging to a file or logging an error.
file_name = os.path.basename(__file__).split(os.path.extsep)[0]
test_logger = logging.getLogger(os.path.extsep.join([__parent_name__,
                                                     file_name]))


class TestGetTopLogger(ut.TestCase):

    def test_logger_class(self):
        self.assertIsInstance(top_logger, RaisableLogger)
        self.assertIs(get_top_logger(), top_logger)


class TestRestoreConfig(ut.TestCase):

    @restore_config
    def test_restore_config_no_error(self):
        def config_modifier_no_restore():
            set_config(verbose=2)
            self.assertEqual(get_config().verbose, 2)
            self.assertEqual(get_config().quiet, 0)

        @restore_config
        def config_modifier_restore():
            set_config(verbose=2)
            self.assertEqual(get_config().verbose, 2)
            self.assertEqual(get_config().quiet, 0)

        set_config(quiet=1)
        config_modifier_restore()
        self.assertEqual(get_config().verbose, 0)
        self.assertEqual(get_config().quiet, 1)
        config_modifier_no_restore()
        self.assertEqual(get_config().verbose, 2)
        self.assertEqual(get_config().quiet, 0)

    @restore_config
    def test_restore_config_error(self):
        def config_modifier_no_restore():
            set_config(verbose=2)
            self.assertEqual(get_config().verbose, 2)
            self.assertEqual(get_config().quiet, 0)
            return 1 / 0

        @restore_config
        def config_modifier_restore():
            set_config(verbose=2)
            self.assertEqual(get_config().verbose, 2)
            self.assertEqual(get_config().quiet, 0)
            return 1 / 0

        set_config(quiet=1)
        self.assertRaises(ZeroDivisionError, config_modifier_restore)
        self.assertEqual(get_config().verbose, 0)
        self.assertEqual(get_config().quiet, 1)
        self.assertRaises(ZeroDivisionError, config_modifier_no_restore)
        self.assertEqual(get_config().verbose, 2)
        self.assertEqual(get_config().quiet, 0)


class TestGetVerbosity(ut.TestCase):

    @restore_config
    def test_get_verbosity(self):
        set_config(quiet=1)
        self.assertEqual(get_verbosity(0, 0), logging.WARNING)
        self.assertEqual(get_verbosity(1, 0), logging.INFO)
        self.assertEqual(get_verbosity(2, 0), logging.DEBUG)
        self.assertEqual(get_verbosity(3, 0), logging.DEBUG)
        self.assertEqual(get_verbosity(0, 1), logging.ERROR)
        self.assertEqual(get_verbosity(0, 2), logging.CRITICAL)
        self.assertEqual(get_verbosity(0, 3), logging.CRITICAL)
        self.assertEqual(get_verbosity(1, 1), logging.WARNING)
        self.assertEqual(get_verbosity(2, 1), logging.WARNING)
        self.assertEqual(get_verbosity(1, 2), logging.WARNING)


class TestExcInfo(ut.TestCase):

    @restore_config
    def test_exc_info(self):
        set_config()
        self.assertFalse(exc_info())
        set_config(verbose=1)
        self.assertFalse(exc_info())
        set_config(verbose=2)
        self.assertTrue(exc_info())
        set_config(quiet=1)
        self.assertFalse(exc_info())
        set_config(quiet=2)
        self.assertFalse(exc_info())


class TestLoggingRaiseOnError(ut.TestCase):

    @restore_config
    def test_raise_on_error(self):
        set_config(raise_on_error=True)
        self.assertRaisesRegex(RuntimeError,
                               "An error has occurred",
                               test_logger.error,
                               "An error has occurred")
        self.assertRaisesRegex(RuntimeError,
                               "A critical error has occurred",
                               test_logger.critical,
                               "A critical error has occurred")
        self.assertRaisesRegex(ZeroDivisionError,
                               "Cannot divide by 0",
                               test_logger.error,
                               ZeroDivisionError("Cannot divide by 0"))
        self.assertRaisesRegex(ZeroDivisionError,
                               "Cannot divide by 0",
                               test_logger.critical,
                               ZeroDivisionError("Cannot divide by 0"))

    @restore_config
    def test_no_raise_on_error(self):
        set_config(raise_on_error=False)
        # Do not log any errors, even critical errors.
        top_logger.setLevel(logging.CRITICAL + 1)
        # None of these calls should raise an error.
        test_logger.error("An error has occurred")
        test_logger.critical("A critical error has occurred")
        test_logger.error(ZeroDivisionError("Cannot divide by 0"))
        test_logger.critical(ZeroDivisionError("Cannot divide by 0"))


class TestEraseConfig(ut.TestCase):

    def test_erase_config(self):
        set_config(verbose=1, raise_on_error=True)
        self.assertEqual(top_logger.level, logging.DEBUG)
        self.assertTrue(top_logger.raise_on_error)
        self.assertEqual(len(top_logger.handlers), 1)
        erase_config()
        self.assertEqual(top_logger.level, logging.WARNING)
        self.assertFalse(top_logger.raise_on_error)
        self.assertEqual(len(top_logger.handlers), 0)


class TestSetConfig(ut.TestCase):

    @restore_config
    def test_defaults(self):
        erase_config()
        self.assertListEqual(top_logger.handlers, [])
        set_config()
        self.assertFalse(top_logger.raise_on_error)
        self.assertEqual(top_logger.level, logging.DEBUG)
        self.assertEqual(len(top_logger.handlers), 1)
        stream_handler = top_logger.handlers[0]
        self.assertIsInstance(stream_handler, logging.StreamHandler)
        self.assertEqual(stream_handler.level, logging.WARNING)
        self.assertIsInstance(stream_handler.formatter, ColorFormatter)

    @restore_config
    def test_verbose(self):
        for verbose, level in {1: logging.INFO, 2: logging.DEBUG}.items():
            set_config(verbose=verbose)
            self.assertEqual(top_logger.level, logging.DEBUG)
            self.assertEqual(len(top_logger.handlers), 1)
            stream_handler = top_logger.handlers[0]
            self.assertIsInstance(stream_handler, logging.StreamHandler)
            self.assertEqual(stream_handler.level, level)

    @restore_config
    def test_quiet(self):
        for quiet, level in {1: logging.ERROR, 2: logging.CRITICAL}.items():
            set_config(quiet=quiet)
            self.assertEqual(top_logger.level, logging.DEBUG)
            self.assertEqual(len(top_logger.handlers), 1)
            stream_handler = top_logger.handlers[0]
            self.assertIsInstance(stream_handler, logging.StreamHandler)
            self.assertEqual(stream_handler.level, level)

    @restore_config
    def test_log_file(self):
        log_file = "test.log"
        msg = "Some logging text"
        try:
            set_config(log_file=log_file)
            self.assertEqual(top_logger.level, logging.DEBUG)
            self.assertEqual(len(top_logger.handlers), 2)
            stream_handler = top_logger.handlers[0]
            file_handler = top_logger.handlers[1]
            self.assertIsInstance(stream_handler, logging.StreamHandler)
            self.assertEqual(stream_handler.level, logging.WARNING)
            self.assertIsInstance(stream_handler.formatter, ColorFormatter)
            self.assertIsInstance(file_handler, logging.FileHandler)
            self.assertEqual(file_handler.level, logging.DEBUG)
            self.assertEqual(file_handler.baseFilename,
                             os.path.abspath(log_file))
            self.assertIsInstance(file_handler.formatter, logging.Formatter)
            self.assertNotIsInstance(file_handler.formatter, ColorFormatter)
            # Test logging a message to the file.
            with open(log_file) as f:
                self.assertEqual(f.readlines(), [])
            test_logger.log(logging.DEBUG, msg)
            with open(log_file) as f:
                lines = f.readlines()
                self.assertEqual(len(lines), 3)
                self.assertTrue(lines[0].startswith("LOGMSG>\t"))
                self.assertEqual(lines[1], f"{msg}\n")
                self.assertEqual(lines[2], "\n")
        finally:
            try:
                os.remove(log_file)
            except FileNotFoundError:
                pass

    @restore_config
    def test_no_log_color(self):
        log_file = "test.log"
        try:
            set_config(log_file=log_file, log_color=False)
            self.assertEqual(top_logger.level, logging.DEBUG)
            self.assertEqual(len(top_logger.handlers), 2)
            stream_handler = top_logger.handlers[0]
            file_handler = top_logger.handlers[1]
            self.assertIsInstance(stream_handler, logging.StreamHandler)
            self.assertEqual(stream_handler.level, logging.WARNING)
            self.assertIsInstance(stream_handler.formatter, logging.Formatter)
            self.assertNotIsInstance(stream_handler.formatter, ColorFormatter)
            self.assertIsInstance(file_handler, logging.FileHandler)
            self.assertEqual(file_handler.level, logging.DEBUG)
            self.assertEqual(file_handler.baseFilename,
                             os.path.abspath(log_file))
            self.assertIsInstance(file_handler.formatter, logging.Formatter)
            self.assertNotIsInstance(file_handler.formatter, ColorFormatter)
        finally:
            try:
                os.remove(log_file)
            except FileNotFoundError:
                pass

    @restore_config
    def test_raise_on_error(self):
        set_config(raise_on_error=True)
        self.assertTrue(top_logger.raise_on_error)


class TestGetConfig(ut.TestCase):

    @restore_config
    def test_get_config(self):
        log_file = os.path.abspath("test.log")
        options = {"verbose": [0, 1, 2],
                   "quiet": [0, 1, 2],
                   "log_file": [None, log_file],
                   "log_color": [True, False],
                   "raise_on_error": [False, True]}
        try:
            for choice in product(*options.values()):
                config = dict(zip(options.keys(), choice, strict=True))
                if config["verbose"] and config["quiet"]:
                    # Cannot set verbose and quiet simultaneously.
                    continue
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
