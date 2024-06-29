import unittest as ut
from pathlib import Path

from seismicrna.core.path import randdir
from seismicrna.relate.sam import (line_attrs,
                                   _iter_records_paired)


class TestLineAttrs(ut.TestCase):

    def test_single(self):
        line = "readname\t0\trefname\t1\t42\t3=\t=\t115\t278\tGAA\tFFF"
        expect = "readname", False, False
        self.assertEqual(line_attrs(line), expect)

    def test_paired_improper(self):
        line = "readname\t1\trefname\t1\t42\t3=\t=\t115\t278\tGAA\tFFF"
        expect = "readname", True, False
        self.assertEqual(line_attrs(line), expect)

    def test_paired_proper(self):
        line = "readname\t3\trefname\t1\t42\t3=\t=\t115\t278\tGAA\tFFF"
        expect = "readname", True, True
        self.assertEqual(line_attrs(line), expect)


def write_sam(text: str):
    sam_file = randdir().joinpath("test.sam")
    with open(sam_file, "w") as f:
        f.write(text)
    return sam_file


def delete_sam(sam_file: Path):
    sam_file.unlink()
    sam_file.parent.rmdir()


class TestIterRecordsPaired(ut.TestCase):

    def run_test_valid(self, lines: list[str], expect: list[tuple[str, str]]):
        text = "\n".join(lines)
        sam_file = write_sam(text)
        try:
            result = [(rec1.rstrip(), rec2.rstrip())
                      for rec1, rec2
                      in _iter_records_paired(sam_file, 0, len(text))]
            self.assertEqual(result, expect)
        except Exception:
            raise
        finally:
            delete_sam(sam_file)

    def run_test_invalid(self, lines: list[str], expect: str):
        text = "\n".join(lines)
        sam_file = write_sam(text)
        try:
            self.assertRaisesRegex(
                ValueError,
                expect,
                lambda: list(_iter_records_paired(sam_file, 0, len(text)))
            )
        except Exception:
            raise
        finally:
            delete_sam(sam_file)

    def test_blank(self):
        lines = []
        expect = []
        self.run_test_valid(lines, expect)

    def test_one_single(self):
        lines = ["read1\t0\ta"]
        expect = "Read 'read1' in .+ is not paired-end"
        self.run_test_invalid(lines, expect)

    def test_one_improper(self):
        lines = ["read1\t1\ta"]
        expect = [("read1\t1\ta", "read1\t1\ta")]
        self.run_test_valid(lines, expect)

    def test_one_proper(self):
        lines = ["read1\t3\ta"]
        expect = ("Read 'read1' in .+ is properly paired but has no mate, "
                  "which indicates a bug")
        self.run_test_invalid(lines, expect)

    def test_two_mated_improper(self):
        lines = ["read1\t1\ta",
                 "read1\t1\tb"]
        expect = [("read1\t1\ta", "read1\t1\ta"),
                  ("read1\t1\tb", "read1\t1\tb")]
        self.run_test_valid(lines, expect)

    def test_two_mated_proper(self):
        lines = ["read1\t3\ta",
                 "read1\t3\tb"]
        expect = [("read1\t3\ta", "read1\t3\tb")]
        self.run_test_valid(lines, expect)

    def test_two_mated_improper_1(self):
        lines = ["read1\t1\ta",
                 "read1\t3\tb"]
        expect = ("Read 'read1' in .+ has only one properly paired mate, "
                  "which indicates a bug")
        self.run_test_invalid(lines, expect)

    def test_two_mated_improper_2(self):
        lines = ["read1\t3\ta",
                 "read1\t1\tb"]
        expect = ("Read 'read1' in .+ has only one properly paired mate, "
                  "which indicates a bug")
        self.run_test_invalid(lines, expect)

    def test_two_unmated_proper(self):
        lines = ["read1\t3\ta",
                 "read2\t3\tb"]
        expect = ("Read 'read1' in .+ is properly paired but has no mate, "
                  "which indicates a bug")
        self.run_test_invalid(lines, expect)

    def test_two_unmated_improper(self):
        lines = ["read1\t1\ta",
                 "read2\t1\tb"]
        expect = [("read1\t1\ta", "read1\t1\ta"),
                  ("read2\t1\tb", "read2\t1\tb")]
        self.run_test_valid(lines, expect)


if __name__ == "__main__":
    ut.main()
