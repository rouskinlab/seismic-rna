import os
import unittest as ut
from itertools import product
from pathlib import Path

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
        self.assertEqual(Level.SEVERE, -3)
        self.assertEqual(Level.ERROR, -2)
        self.assertEqual(Level.WARNING, -1)
        self.assertEqual(Level.COMMAND, 0)
        self.assertEqual(Level.TASK, 1)
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
        for verbosity in [Level.COMMAND,
                          Level.WARNING,
                          Level.ERROR,
                          Level.SEVERE]:
            set_config(verbosity=verbosity)
            self.assertFalse(exc_info())
        for verbosity in [Level.DETAIL,
                          Level.ROUTINE,
                          Level.TASK]:
            set_config(verbosity=verbosity)
            self.assertTrue(exc_info())


class TestLoggingRaiseOnError(ut.TestCase):

    @restore_config
    def test_raise_on_error(self):
        set_config(verbosity=(Level.SEVERE - 1), raise_on_error=True)
        logger.warning("An error has occurred")
        self.assertRaisesRegex(RuntimeError,
                               "An error has occurred",
                               logger.error,
                               "An error has occurred")
        self.assertRaisesRegex(RuntimeError,
                               "A fatal error has occurred",
                               logger.severe,
                               "A fatal error has occurred")
        logger.warning(ZeroDivisionError("Cannot divide by 0"))
        self.assertRaisesRegex(ZeroDivisionError,
                               "Cannot divide by 0",
                               logger.error,
                               ZeroDivisionError("Cannot divide by 0"))
        self.assertRaisesRegex(ZeroDivisionError,
                               "Cannot divide by 0",
                               logger.severe,
                               ZeroDivisionError("Cannot divide by 0"))

    @restore_config
    def test_no_raise_on_error(self):
        set_config(verbosity=(Level.SEVERE - 1), raise_on_error=False)
        # None of these calls should raise an error.
        logger.error("An error has occurred")
        logger.severe("A fatal error has occurred")
        logger.error(ZeroDivisionError("Cannot divide by 0"))
        logger.severe(ZeroDivisionError("Cannot divide by 0"))


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
        log_file_path = Path(os.path.abspath("test.log"))
        msg1 = "Some logging text"
        msg2 = "More logging text"
        try:
            set_config(log_file_path=log_file_path)
            self.assertEqual(logger.file_stream.filterer.verbosity,
                             FILE_VERBOSITY)
            self.assertIs(logger.file_stream.formatter.formatter,
                          format_logfile)
            self.assertEqual(logger.file_stream.file_path, log_file_path)
            # Test logging a message to the file.
            logger.routine(msg1)
            logger.detail(msg2)
            # The file stream must be closed explicitly to flush the log
            # messages to the file before reading with readlines().
            logger.file_stream.close()
            with open(log_file_path) as log_file:
                lines = log_file.readlines()
                self.assertEqual(len(lines), 6)
                self.assertTrue(
                    lines[0].startswith(f"LOGMSG> {Level.ROUTINE.name}")
                )
                self.assertEqual(lines[1], f"{msg1}\n")
                self.assertEqual(lines[2], "\n")
                self.assertTrue(
                    lines[3].startswith(f"LOGMSG> {Level.DETAIL.name}")
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
        log_file_path = Path(os.path.abspath("test.log"))
        options = {"verbosity": list(Level),
                   "log_file_path": [None, log_file_path],
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
                os.remove(log_file_path)
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    ut.main()

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
