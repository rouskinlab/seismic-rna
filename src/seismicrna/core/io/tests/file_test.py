import unittest as ut
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp, mkstemp

import numpy as np

from seismicrna.core.io.file import (make_temp_backup,
                                     restore_temp_backup)

rng = np.random.default_rng()

FILE_SIZE = 64
FILE_SUFFIX = ".file"


class SourceBackupDirs(object):

    def __init__(self):
        self._source_dir = None
        self._backup_dir = None

    @property
    def source_dir(self) -> Path:
        if self._source_dir is None:
            raise TypeError("Source directory is not set")
        return self._source_dir

    @property
    def backup_dir(self) -> Path:
        if self._backup_dir is None:
            raise TypeError("Backup directory is not set")
        return self._backup_dir

    def __enter__(self):
        # Make temporary source and backup directories.
        self._source_dir = Path(mkdtemp()).resolve(strict=True)
        self._backup_dir = Path(mkdtemp()).resolve(strict=True)
        return self.source_dir, self.backup_dir

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Delete the temporary source and backup directories.
        rmtree(self.source_dir, ignore_errors=True)
        rmtree(self.backup_dir, ignore_errors=True)
        self._source_dir = None
        self._backup_dir = None


class TestMakeTempBackup(ut.TestCase):

    @staticmethod
    def make_file(in_dir: Path):
        _, path = mkstemp(dir=in_dir, suffix=FILE_SUFFIX)
        with open(path, "wb") as f:
            f.write(content := rng.bytes(FILE_SIZE))
        return Path(path), content

    @staticmethod
    def make_dir(in_dir: Path):
        return Path(mkdtemp(dir=in_dir))

    def test_backup_file(self):
        with SourceBackupDirs() as (source_dir, backup_dir):
            source_file, content = self.make_file(source_dir)
            backup_file = make_temp_backup(source_file,
                                           source_dir,
                                           backup_dir)
            self.assertEqual(backup_file.parent, backup_dir)
            self.assertEqual(backup_file.name, source_file.name)
            with open(backup_file, "rb") as f:
                self.assertEqual(f.read(), content)

    def test_backup_dir(self):
        with SourceBackupDirs() as (source_dir, backup_dir):
            source_subdir = self.make_dir(source_dir)
            source_file1, content1 = self.make_file(source_subdir)
            source_file2, content2 = self.make_file(source_subdir)
            backup_subdir = make_temp_backup(source_subdir,
                                             source_dir,
                                             backup_dir)
            self.assertEqual(backup_subdir,
                             backup_dir.joinpath(source_subdir.name))
            backup_file1 = backup_subdir.joinpath(source_file1.name)
            backup_file2 = backup_subdir.joinpath(source_file2.name)
            with open(backup_file1, "rb") as f:
                self.assertEqual(f.read(), content1)
            with open(backup_file2, "rb") as f:
                self.assertEqual(f.read(), content2)

    def test_restore_file(self):
        with SourceBackupDirs() as (source_dir, backup_dir):
            source_file, content = self.make_file(source_dir)
            backup_file = make_temp_backup(source_file,
                                           source_dir,
                                           backup_dir)
            self.assertTrue(source_file.is_file())
            source_file.unlink()
            self.assertFalse(source_file.is_file())
            self.assertEqual(restore_temp_backup(source_file,
                                                 source_dir,
                                                 backup_dir),
                             backup_file)
            self.assertTrue(source_file.is_file())
            with open(source_file, "rb") as f:
                self.assertEqual(f.read(), content)

    def test_restore_dir(self):
        with SourceBackupDirs() as (source_dir, backup_dir):
            source_subdir = self.make_dir(source_dir)
            source_file1, content1 = self.make_file(source_subdir)
            source_file2, content2 = self.make_file(source_subdir)
            backup_subdir = make_temp_backup(source_subdir,
                                             source_dir,
                                             backup_dir)
            self.assertTrue(source_subdir.is_dir())
            self.assertTrue(source_file1.is_file())
            self.assertTrue(source_file2.is_file())
            rmtree(source_subdir)
            self.assertFalse(source_subdir.is_dir())
            self.assertFalse(source_file1.is_file())
            self.assertFalse(source_file2.is_file())
            self.assertEqual(restore_temp_backup(source_subdir,
                                                 source_dir,
                                                 backup_dir),
                             backup_subdir)
            self.assertTrue(source_subdir.is_dir())
            self.assertTrue(source_file1.is_file())
            self.assertTrue(source_file2.is_file())
            with open(source_file1, "rb") as f:
                self.assertEqual(f.read(), content1)
            with open(source_file2, "rb") as f:
                self.assertEqual(f.read(), content2)


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
