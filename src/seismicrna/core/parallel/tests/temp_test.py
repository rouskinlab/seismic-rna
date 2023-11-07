import unittest as ut
from logging import Filter, LogRecord
from pathlib import Path
from tempfile import mkdtemp

from ..temp import LOCK_DIR, lock_temp_dir, logger as temp_logger


def get_lock(temp_dir: Path):
    return temp_dir.joinpath(LOCK_DIR)


def lock(temp_dir: Path):
    lock_dir = get_lock(temp_dir)
    lock_dir.mkdir(parents=False, exist_ok=False)
    return lock_dir


def unlock(temp_dir: Path):
    try:
        get_lock(temp_dir).rmdir()
    except FileNotFoundError:
        pass


def make_temp(*args, **kwargs):
    """ Make a new temporary directory with a random path. """
    return Path(mkdtemp(*args, **kwargs))


def rm_temp(temp_dir: Path):
    """ Remove a temporary directory, if it exists. """
    unlock(temp_dir)
    try:
        temp_dir.rmdir()
    except FileNotFoundError:
        return False
    return True


def name_temp(*args, **kwargs):
    """ Get the path of a temporary directory that does not exist. """
    temp_dir = make_temp(*args, **kwargs)
    rm_temp(temp_dir)
    return temp_dir


def make_lock_temp(*args, **kwargs):
    """ Make and lock a new temporary directory. """
    temp_dir = make_temp(*args, **kwargs)
    return temp_dir, lock(temp_dir)


def name_lock_temp(*args, **kwargs):
    temp_dir = name_temp(*args, **kwargs)
    return temp_dir, get_lock(temp_dir)


@lock_temp_dir
def run_func(*_, temp_dir: Path, keep_temp: bool, **__):
    """ Placeholder for run() function. """
    return keep_temp, temp_dir.is_dir(), get_lock(temp_dir).is_dir()


class TestLockTempDir(ut.TestCase):
    """ Test decorator `lock_temp_dir`. """

    class LockErrFilter(Filter):

        def filter(self, record: LogRecord):
            return "currently being used" not in record.msg

    class TempErrFilter(Filter):

        def filter(self, record: LogRecord):
            return "If any needed files" not in record.msg

    def test_wraps(self):
        self.assertEqual(run_func.__name__, "run_func")
        self.assertEqual(run_func.__doc__, " Placeholder for run() function. ")

    def test_new_keep_temp(self):
        temp_dir, lock_dir = name_lock_temp()
        try:
            # The directory should not exist initially.
            self.assertFalse(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The directory should be created by lock_temp_dir().
            self.assertEqual(run_func(temp_dir=temp_dir, keep_temp=True),
                             (True, True, True))
            # The directory should not be deleted by lock_temp_dir().
            self.assertTrue(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            rm_temp(temp_dir)

    def test_new_erase_temp(self):
        temp_dir, lock_dir = name_lock_temp()
        try:
            # The directory should not exist initially.
            self.assertFalse(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The directory should be created by lock_temp_dir().
            self.assertEqual(run_func(temp_dir=temp_dir, keep_temp=False),
                             (False, True, True))
            # The directory should be deleted by lock_temp_dir().
            self.assertFalse(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            rm_temp(temp_dir)

    def test_exists_keep_temp(self):
        temp_dir = make_temp()
        lock_dir = get_lock(temp_dir)
        try:
            # The directory should exist initially.
            self.assertTrue(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The function should run normally.
            self.assertEqual(run_func(temp_dir=temp_dir, keep_temp=True),
                             (True, True, True))
            # The directory should not be deleted by lock_temp_dir().
            self.assertTrue(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            rm_temp(temp_dir)

    def test_exists_erase_temp(self):
        temp_dir = make_temp()
        lock_dir = get_lock(temp_dir)
        temp_logger.addFilter(temp_err := self.TempErrFilter())
        try:
            # The directory should exist initially.
            self.assertTrue(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The function should fail.
            self.assertRaises(SystemExit, run_func,
                              temp_dir=temp_dir, keep_temp=False)
            # The directory should not be deleted by lock_temp_dir.
            self.assertTrue(temp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            temp_logger.removeFilter(temp_err)
            rm_temp(temp_dir)

    def test_locked(self):
        temp_dir, lock_dir = make_lock_temp()
        temp_logger.addFilter(lock_err := self.LockErrFilter())
        try:
            for keep_temp in (True, False):
                # The directory should exist initially.
                self.assertTrue(temp_dir.is_dir())
                self.assertTrue(lock_dir.is_dir())
                # The function should fail.
                self.assertRaises(SystemExit, run_func,
                                  temp_dir=temp_dir, keep_temp=keep_temp)
                # The directory should not be deleted by lock_temp_dir.
                self.assertTrue(temp_dir.is_dir())
                self.assertTrue(lock_dir.is_dir())
        finally:
            temp_logger.removeFilter(lock_err)
            rm_temp(temp_dir)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
