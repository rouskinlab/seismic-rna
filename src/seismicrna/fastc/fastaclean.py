"""

FASTA Cleaner Module

"""

import re
from logging import getLogger
from pathlib import Path

from ..core.fasta import extract_fasta_seqname, format_fasta_name_line
from ..core.seq import Seq, BASEN

logger = getLogger(__name__)


def get_non_seq_regex(seq_type: type[Seq]):
    return re.compile("[" + "".join(seq_type.get_nonalphaset() + {'\n'}) + "]")


class FastaCleaner(object):
    __slots__ = "non_seq_regex",

    def __init__(self, seq_type: type[Seq]):
        self.non_seq_regex = get_non_seq_regex(seq_type)

    def run(self, ifasta: Path, ofasta: Path, force: bool = False):
        with open(ifasta) as fi, open(ofasta, 'w' if force else 'x') as fo:
            for line in fi:
                fo.write(format_fasta_name_line(name)
                         if (name := extract_fasta_seqname(line))
                         else self.non_seq_regex.sub(BASEN, line.upper()))
