import os
import tempfile
import unittest as ut

import click
from click.testing import CliRunner

from seismicrna.core.arg.glob import GlobPath, expand_ct_pos_5, flatten_glob_results


def _make_multi_cli(**kwargs):
    @click.command()
    @click.option(
        "-x",
        type=GlobPath(exists=True),
        multiple=True,
        default=(),
        callback=flatten_glob_results,
        **kwargs,
    )
    def cmd(x):
        for path in x:
            click.echo(path)

    return cmd


def _make_ct_pos_5_cli():
    @click.command()
    @click.option(
        "--ct-pos-5",
        type=(GlobPath(exists=True), int),
        multiple=True,
        default=(),
        callback=expand_ct_pos_5,
    )
    def cmd(ct_pos_5):
        for path, pos in ct_pos_5:
            click.echo(f"{path}\t{pos}")

    return cmd


class TestGlobPathMultiple(ut.TestCase):
    def test_literal_path_passes_through(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = os.path.join(tmp, "a.fastq")
            open(f, "w").close()
            runner = CliRunner()
            result = runner.invoke(_make_multi_cli(), ["-x", f])
            self.assertEqual(result.exit_code, 0)
            self.assertEqual(result.output.strip().splitlines(), [f])

    def test_glob_expands_multiple_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            f1 = os.path.join(tmp, "s_R1.fastq")
            f2 = os.path.join(tmp, "s_R2.fastq")
            open(f1, "w").close()
            open(f2, "w").close()
            pattern = os.path.join(tmp, "s_R[12].fastq")
            runner = CliRunner()
            result = runner.invoke(_make_multi_cli(), ["-x", pattern])
            self.assertEqual(result.exit_code, 0)
            self.assertEqual(result.output.strip().splitlines(), sorted([f1, f2]))

    def test_glob_combined_with_literal(self):
        with tempfile.TemporaryDirectory() as tmp:
            f1 = os.path.join(tmp, "a.fastq")
            f2 = os.path.join(tmp, "b.fastq")
            f3 = os.path.join(tmp, "c.fastq")
            for f in (f1, f2, f3):
                open(f, "w").close()
            pattern = os.path.join(tmp, "[ab].fastq")
            runner = CliRunner()
            result = runner.invoke(_make_multi_cli(), ["-x", pattern, "-x", f3])
            self.assertEqual(result.exit_code, 0)
            self.assertEqual(result.output.strip().splitlines(), [f1, f2, f3])

    def test_no_match_raises_file_not_found(self):
        with tempfile.TemporaryDirectory() as tmp:
            pattern = os.path.join(tmp, "nope_*.fastq")
            runner = CliRunner()
            result = runner.invoke(_make_multi_cli(), ["-x", pattern])
            self.assertNotEqual(result.exit_code, 0)
            self.assertIsInstance(result.exception, FileNotFoundError)
            self.assertIn("No files matched", str(result.exception))


class TestCtPos5Expansion(ut.TestCase):
    def test_pattern_paired_with_position(self):
        with tempfile.TemporaryDirectory() as tmp:
            f1 = os.path.join(tmp, "a.ct")
            f2 = os.path.join(tmp, "b.ct")
            open(f1, "w").close()
            open(f2, "w").close()
            pattern = os.path.join(tmp, "*.ct")
            runner = CliRunner()
            result = runner.invoke(_make_ct_pos_5_cli(), ["--ct-pos-5", pattern, "7"])
            self.assertEqual(result.exit_code, 0)
            self.assertEqual(
                sorted(result.output.strip().splitlines()), [f"{f1}\t7", f"{f2}\t7"]
            )

    def test_literal_passes_through(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = os.path.join(tmp, "a.ct")
            open(f, "w").close()
            runner = CliRunner()
            result = runner.invoke(_make_ct_pos_5_cli(), ["--ct-pos-5", f, "3"])
            self.assertEqual(result.exit_code, 0)
            self.assertEqual(result.output.strip(), f"{f}\t3")


if __name__ == "__main__":
    ut.main(verbosity=2)
