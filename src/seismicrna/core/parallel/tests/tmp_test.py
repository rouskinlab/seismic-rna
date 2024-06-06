import unittest as ut
from logging import Filter, LogRecord
from pathlib import Path
from tempfile import mkdtemp

from seismicrna.core.parallel.tmp import (LOCK_DIR,
                                          lock_tmp_dir,
                                          logger as tmp_logger)


def get_lock(tmp_dir: Path):
    return tmp_dir.joinpath(LOCK_DIR)


def lock(tmp_dir: Path):
    lock_dir = get_lock(tmp_dir)
    lock_dir.mkdir(parents=False, exist_ok=False)
    return lock_dir


def unlock(tmp_dir: Path):
    try:
        get_lock(tmp_dir).rmdir()
    except FileNotFoundError:
        pass


def make_tmp(*args, **kwargs):
    """ Make a new temporary directory with a random path. """
    return Path(mkdtemp(*args, **kwargs))


def rm_tmp(tmp_dir: Path):
    """ Remove a temporary directory, if it exists. """
    unlock(tmp_dir)
    try:
        tmp_dir.rmdir()
    except FileNotFoundError:
        return False
    return True


def name_tmp(*args, **kwargs):
    """ Get the path of a temporary directory that does not exist. """
    tmp_dir = make_tmp(*args, **kwargs)
    rm_tmp(tmp_dir)
    return tmp_dir


def make_lock_tmp(*args, **kwargs):
    """ Make and lock a new temporary directory. """
    tmp_dir = make_tmp(*args, **kwargs)
    return tmp_dir, lock(tmp_dir)


def name_lock_tmp(*args, **kwargs):
    tmp_dir = name_tmp(*args, **kwargs)
    return tmp_dir, get_lock(tmp_dir)


@lock_tmp_dir
def run_func(*_, tmp_dir: Path, keep_tmp: bool, **__):
    """ Placeholder for run() function. """
    return keep_tmp, tmp_dir.is_dir(), get_lock(tmp_dir).is_dir()


class TestLockTmpDir(ut.TestCase):
    """ Test decorator `lock_tmp_dir`. """

    class LockErrFilter(Filter):

        def filter(self, record: LogRecord):
            return "currently being used" not in record.msg

    class TmpErrFilter(Filter):

        def filter(self, record: LogRecord):
            return "Please either delete it" not in record.msg

    def test_wraps(self):
        self.assertEqual(run_func.__name__, "run_func")
        self.assertEqual(run_func.__doc__, " Placeholder for run() function. ")

    def test_new_keep_tmp(self):
        tmp_dir, lock_dir = name_lock_tmp()
        try:
            # The directory should not exist initially.
            self.assertFalse(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The directory should be created by lock_tmp_dir().
            self.assertEqual(run_func(tmp_dir=tmp_dir, keep_tmp=True),
                             (True, True, True))
            # The directory should not be deleted by lock_tmp_dir().
            self.assertTrue(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            rm_tmp(tmp_dir)

    def test_new_erase_tmp(self):
        tmp_dir, lock_dir = name_lock_tmp()
        try:
            # The directory should not exist initially.
            self.assertFalse(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The directory should be created by lock_tmp_dir().
            self.assertEqual(run_func(tmp_dir=tmp_dir, keep_tmp=False),
                             (False, True, True))
            # The directory should be deleted by lock_tmp_dir().
            self.assertFalse(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            rm_tmp(tmp_dir)

    def test_exists_keep_tmp(self):
        tmp_dir = make_tmp()
        lock_dir = get_lock(tmp_dir)
        try:
            # The directory should exist initially.
            self.assertTrue(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The function should run normally.
            self.assertEqual(run_func(tmp_dir=tmp_dir, keep_tmp=True),
                             (True, True, True))
            # The directory should not be deleted by lock_tmp_dir().
            self.assertTrue(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            rm_tmp(tmp_dir)

    def test_exists_erase_tmp(self):
        tmp_dir = make_tmp()
        lock_dir = get_lock(tmp_dir)
        tmp_logger.addFilter(tmp_err := self.TmpErrFilter())
        try:
            # The directory should exist initially.
            self.assertTrue(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
            # The function should fail.
            self.assertRaises(SystemExit, run_func,
                              tmp_dir=tmp_dir, keep_tmp=False)
            # The directory should not be deleted by lock_tmp_dir.
            self.assertTrue(tmp_dir.is_dir())
            self.assertFalse(lock_dir.is_dir())
        finally:
            tmp_logger.removeFilter(tmp_err)
            rm_tmp(tmp_dir)

    def test_locked(self):
        tmp_dir, lock_dir = make_lock_tmp()
        tmp_logger.addFilter(lock_err := self.LockErrFilter())
        try:
            for keep_tmp in (True, False):
                # The directory should exist initially.
                self.assertTrue(tmp_dir.is_dir())
                self.assertTrue(lock_dir.is_dir())
                # The function should fail.
                self.assertRaises(SystemExit, run_func,
                                  tmp_dir=tmp_dir, keep_tmp=keep_tmp)
                # The directory should not be deleted by lock_tmp_dir.
                self.assertTrue(tmp_dir.is_dir())
                self.assertTrue(lock_dir.is_dir())
        finally:
            tmp_logger.removeFilter(lock_err)
            rm_tmp(tmp_dir)


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
