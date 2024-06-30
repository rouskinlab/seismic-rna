import unittest as ut
from itertools import product
from os import linesep, remove
from pathlib import Path
from string import printable, whitespace
from tempfile import NamedTemporaryFile as NTFile

from seismicrna.core import path
from seismicrna.core.seq.fasta import (FASTA_NAME_MARK,
                                       FASTA_NAME_CHARS,
                                       valid_fasta_seqname,
                                       format_fasta_name_line,
                                       format_fasta_record,
                                       format_fasta_seq_lines,
                                       parse_fasta,
                                       write_fasta)
from seismicrna.core.seq.xna import DNA, RNA


class TestValidFastaSeqname(ut.TestCase):
    """ Test valid_fasta_seqname. """

    def test_name_mark(self):
        self.assertEqual(FASTA_NAME_MARK, '>')

    def test_valid(self):
        """ Test that it works on valid name lines. """
        for name in FASTA_NAME_CHARS:
            prefix = f"{FASTA_NAME_MARK}{name}"
            for line in [prefix, f"{prefix}\n"]:
                self.assertEqual(valid_fasta_seqname(line), name)

    def test_misformatted(self):
        """ Test that it fails on misformatted lines. """
        for a, b in product(printable, repeat=2):
            if a != FASTA_NAME_MARK:
                prefix = f"{a}{b}"
                for line in [prefix, f"{prefix}\n"]:
                    self.assertRaisesRegex(ValueError, "is misformatted",
                                           valid_fasta_seqname, line)

    def test_illegal_prefix(self):
        """ Test it fails on names starting with illegal characters. """
        prefixes = set(printable) - set(FASTA_NAME_CHARS)
        for a, b in product(prefixes, FASTA_NAME_CHARS):
            prefix = f"{FASTA_NAME_MARK}{a}{b}"
            for line in [prefix, f"{prefix}\n"]:
                self.assertRaisesRegex(ValueError, "has a blank name",
                                       valid_fasta_seqname, line)

    def test_illegal_suffix(self):
        """ Test it fails on names ending with illegal characters except
        for trailing whitespace, which is simply ignored. """
        suffixes = set(printable) - (set(FASTA_NAME_CHARS) | set(whitespace))
        for a, b in product(FASTA_NAME_CHARS, suffixes):
            prefix = f"{FASTA_NAME_MARK}{a}{b}"
            for line in [prefix, f"{prefix}\n"]:
                self.assertRaisesRegex(ValueError, "has illegal characters",
                                       valid_fasta_seqname, line)

    def test_blank(self):
        """ Test that it fails with blank names. """
        prefixes = [FASTA_NAME_MARK] + [f"{FASTA_NAME_MARK}{w}"
                                        for w in whitespace]
        for prefix in prefixes:
            for line in [prefix, f"{prefix}\n"]:
                self.assertRaisesRegex(ValueError, "has a blank name",
                                       valid_fasta_seqname, line)


class TestFormat(ut.TestCase):
    """ Test format_fasta_name_line and format_fasta_record. """

    def test_format_fasta_name_line(self):
        lines = {
            "name": ">name\n",
            "name ": ">name\n",
            "name\n": ">name\n",
            " name": "> name\n",
        }
        for name, line in lines.items():
            self.assertEqual(format_fasta_name_line(name), line)

    def test_format_fasta_seq_lines(self):
        lines = {
            (DNA("AGTC"), 0): "AGTC\n",
            (DNA("AGTC"), 1): "A\nG\nT\nC\n",
            (DNA("AGTC"), 2): "AG\nTC\n",
            (DNA("AGTC"), 3): "AGT\nC\n",
            (DNA("AGTC"), 4): "AGTC\n",
            (DNA("AGTC"), 5): "AGTC\n",
        }
        for (seq, wrap), line in lines.items():
            self.assertEqual(format_fasta_seq_lines(seq, wrap), line)

    def test_format_fasta_record(self):
        lines = {
            ("name", DNA("ACGTN"), 0): ">name\nACGTN\n",
            ("name", RNA("ACGUN"), 0): ">name\nACGUN\n",
            ("name", DNA("ANCGTN"), 3): ">name\nANC\nGTN\n",
            ("name", RNA("ANCGUN"), 2): ">name\nAN\nCG\nUN\n",
        }
        for (name, seq, wrap), line in lines.items():
            self.assertEqual(format_fasta_record(name, seq, wrap), line)


class TestParseFasta(ut.TestCase):
    """ Test parse_fasta. """

    def test_valid_names(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([">Seq1", "GTACGTGNTCATC",
                                  ">Seq2 ", "AGCTGTGNNT", "ATCG"]))
        try:
            records = list(parse_fasta(filepath, None))
        finally:
            remove(filepath)
        self.assertEqual(records, ["Seq1", "Seq2"])

    def test_valid_dna(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([">Seq1", "GTACGTGNTCATC",
                                  ">Seq2 ", "AGCTGTGNNT", "ATCG"]))
        try:
            records = dict(parse_fasta(filepath, DNA))
        finally:
            remove(filepath)
        self.assertEqual(records, {"Seq1": DNA("GTACGTGNTCATC"),
                                   "Seq2": DNA("AGCTGTGNNTATCG")})

    def test_valid_rna(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([">Seq1", "GUACGUGNUCAUC",
                                  ">Seq2 ", "AGCUGUGNNU", "AUCG"]))
        try:
            records = dict(parse_fasta(filepath, RNA))
        finally:
            remove(filepath)
        self.assertEqual(records, {"Seq1": RNA("GUACGUGNUCAUC"),
                                   "Seq2": RNA("AGCUGUGNNUAUCG")})

    def test_valid_empty(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([]))
        try:
            records = dict(parse_fasta(filepath, DNA))
        finally:
            remove(filepath)
        self.assertEqual(records, {})

    def test_valid_blank_line(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([">Seq1", "GTACGTGNTCATC", linesep, linesep,
                                  ">Seq2 ", linesep, "AGCTGTGNNT", linesep,
                                  "ATCG"]))
        try:
            records = dict(parse_fasta(filepath, DNA))
        finally:
            remove(filepath)
        self.assertEqual(records, {"Seq1": DNA("GTACGTGNTCATC"),
                                   "Seq2": DNA("AGCTGTGNNTATCG")})

    def test_duplicate_name(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([">Seq1", "GTACGTGNTCATC",
                                  ">Seq1 ", "AGCTGTGNNT", "ATCG"]))
        try:
            self.assertRaisesRegex(ValueError,
                                   "Duplicate name 'Seq1'",
                                   lambda: dict(parse_fasta(filepath, DNA)))
        finally:
            remove(filepath)

    def test_invalid_name(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([">Seq1", "GTACGTGNTCATC",
                                  ">Seq2|", "AGCTGTGNNT", "ATCG"]))
        try:
            self.assertRaisesRegex(
                ValueError,
                "FASTA name line '>Seq2|' has illegal characters",
                lambda: dict(parse_fasta(filepath, DNA))
            )
        finally:
            remove(filepath)

    def test_invalid_seq(self):
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=False) as f:
            filepath = Path(f.file.name)
            f.write(linesep.join([">Seq1", "GTACGTGNTCATC ",
                                  ">Seq2", "AGCTGTGNNT", "ATCG"]))
        try:
            self.assertRaisesRegex(ValueError,
                                   "Invalid DNA bases: {' '}",
                                   lambda: dict(parse_fasta(filepath, DNA)))
        finally:
            remove(filepath)


class TestWriteFasta(ut.TestCase):
    """ Test write_fasta. """

    def test_valid_names(self):
        seqs = [("Seq1", DNA("GTACGTGNTCATC")),
                ("Seq2", DNA("AGCTGTGNNTATCG"))]
        text = ">Seq1\nGTACGTGNTCATC\n>Seq2\nAGCTGTGNNTATCG\n"
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=True) as f:
            filepath = Path(f.file.name)
        try:
            write_fasta(filepath, seqs)
            with open(filepath) as f:
                self.assertEqual(f.read(), text)
        finally:
            remove(filepath)

    def test_force(self):
        seqs = [("Seq1", DNA("GTACGTGNTCATC")),
                ("Seq2", DNA("AGCTGTGNNTATCG"))]
        fasta_text = ">Seq1\nGTACGTGNTCATC\n>Seq2\nAGCTGTGNNTATCG\n"
        other_text = "some other text"
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=True) as f:
            filepath = Path(f.file.name)
        try:
            with open(filepath, "w") as f:
                f.write(other_text)
            write_fasta(filepath, seqs, force=True)
            with open(filepath) as f:
                self.assertEqual(f.read(), fasta_text)
        finally:
            remove(filepath)

    def test_blank_name(self):
        seqs = [("", DNA("GTACGTGNTCATC"))]
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=True) as f:
            filepath = Path(f.file.name)
        try:
            self.assertRaisesRegex(ValueError,
                                   "Got blank reference name",
                                   write_fasta,
                                   filepath, seqs)
        finally:
            filepath.unlink(missing_ok=True)

    def test_invalid_name(self):
        seqs = [(" Seq1", DNA("GACGTACTGTACGT"))]
        with NTFile("w", suffix=path.FASTA_EXTS[0], delete=True) as f:
            filepath = Path(f.file.name)
        try:
            self.assertRaisesRegex(
                ValueError,
                "Reference name ' Seq1' has illegal characters",
                write_fasta,
                filepath, seqs
            )
        finally:
            filepath.unlink(missing_ok=True)

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
