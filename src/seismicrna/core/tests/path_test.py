import unittest as ut
from shutil import rmtree

from seismicrna.core.path import (get_seismicrna_source_dir,
                                  randdir,
                                  sanitize,
                                  symlink_if_needed)


class TestGetSeismicRNASourceDir(ut.TestCase):

    def test_get_seismicrna_source_dir(self):
        seismicrna_source_dir = get_seismicrna_source_dir()
        self.assertEqual(seismicrna_source_dir,
                         sanitize(__file__).parent.parent.parent)


class TestSymlinkIfNeeded(ut.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._tmp_dir = None

    def setUp(self):
        self._tmp_dir = randdir()

    def tearDown(self):
        rmtree(self._tmp_dir)
        self._tmp_dir = None

    def test_target_not_exist(self):
        link = self._tmp_dir.joinpath("link")
        target = self._tmp_dir.joinpath("target")
        for link_exists in range(2):
            self.assertRaisesRegex(FileNotFoundError,
                                   str(target),
                                   symlink_if_needed,
                                   link,
                                   target)
            link.mkdir(exist_ok=link_exists)
            self.assertTrue(link.is_dir())

    def test_link_valid(self):
        link = self._tmp_dir.joinpath("link")
        target = self._tmp_dir.joinpath("target")
        target.mkdir()
        for _ in range(2):
            symlink_if_needed(link, target)
            self.assertTrue(link.is_symlink())
            self.assertTrue(link.readlink() == target)

    def test_link_not_symlink(self):
        link = self._tmp_dir.joinpath("link")
        link.mkdir()
        target = self._tmp_dir.joinpath("target")
        target.mkdir()
        self.assertRaisesRegex(OSError,
                               f"{link} is not a symbolic link",
                               symlink_if_needed,
                               link,
                               target)

    def test_link_wrong_symlink(self):
        link = self._tmp_dir.joinpath("link")
        target = self._tmp_dir.joinpath("target")
        target.mkdir()
        target2 = self._tmp_dir.joinpath("target2")
        link.symlink_to(target2)
        self.assertRaisesRegex(
            OSError,
            f"{link} is a symbolic link to {target2}, not to {target}",
            symlink_if_needed,
            link,
            target
        )


if __name__ == "__main__":
    ut.main(verbosity=2)
