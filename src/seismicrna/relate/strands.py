from pathlib import Path

from ..align.write import format_ref_reverse
from ..core.seq import DNA, parse_fasta, write_fasta


def generate_both_strands(ref: str, seq: DNA, rev_label: str):
    """ Yield both the forward and reverse strand for each sequence. """
    yield ref, seq
    yield format_ref_reverse(ref, rev_label), seq.rc


def write_both_strands(fasta_in: Path, fasta_out: Path, rev_label: str):
    """ Write a FASTA file of both forward and reverse strands. """
    write_fasta(fasta_out,
                (strand
                 for ref, seq in parse_fasta(fasta_in, DNA)
                 for strand in generate_both_strands(ref, seq, rev_label)))

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
