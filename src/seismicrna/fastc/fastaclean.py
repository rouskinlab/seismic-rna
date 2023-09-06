"""

FASTA Cleaner Module

"""

import re
from logging import getLogger
from pathlib import Path
from string import printable

from ..core.seq import (Seq, BASEN, FASTA_NAME_REGEX,
                        match_to_seq_name, format_fasta_name_line)

logger = getLogger(__name__)


def get_non_seq_characters(seq_type: type[Seq]):
    return "".join(sorted(set(printable) - (seq_type.get_alphaset() | {'\n'})))


def get_non_seq_regex(seq_type: type[Seq]):
    return re.compile(f"[{get_non_seq_characters(seq_type)}]")


class FastaCleaner(object):
    __slots__ = "non_seq_regex",

    def __init__(self, seq_type: type[Seq]):
        self.non_seq_regex = get_non_seq_regex(seq_type)

    def run(self, ifasta: Path, ofasta: Path, force: bool = False):
        with open(ifasta) as fi, open(ofasta, 'w' if force else 'x') as fo:
            for line in fi:
                fo.write(self.non_seq_regex.sub(BASEN, line.upper())
                         if (match := FASTA_NAME_REGEX.match(line)) is None
                         else format_fasta_name_line(match_to_seq_name(match)))
