"""

FASTA Cleaner Module

"""

import re
from logging import getLogger
from pathlib import Path

from ..core.seq import BASEN, XNA, extract_fasta_seqname, format_fasta_name_line

logger = getLogger(__name__)


def get_non_seq_regex(seq_type: type[XNA]):
    return re.compile("[" + "".join(seq_type.get_nonalphaset() + {'\n'}) + "]")


class FastaCleaner(object):
    __slots__ = "non_seq_regex",

    def __init__(self, seq_type: type[XNA]):
        self.non_seq_regex = get_non_seq_regex(seq_type)

    def run(self, ifasta: Path, ofasta: Path, force: bool = False):
        if force or not ofasta.is_file():
            with open(ifasta) as fi, open(ofasta, 'w' if force else 'x') as fo:
                for line in fi:
                    fo.write(format_fasta_name_line(name)
                             if (name := extract_fasta_seqname(line))
                             else self.non_seq_regex.sub(BASEN, line.upper()))

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
