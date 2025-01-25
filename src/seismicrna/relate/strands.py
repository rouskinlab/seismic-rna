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
