"""
strandedness.py

Author: Matthew Allan
Date: 2020-06-11

Identify whether each alignment in a SAM file came from the plus or
minus strand of the RNA.

Usage:
    - Count the number of plus and minus strands:
    python strandedness.py alignment_map.sam discordant[yes/no]
    - Count and split the plus and minus strands into two new SAM files:
      python strandedness.py alignment_map.sam discordant[yes/no]
      new_file_for_plus_alignments.sam new_file_for_minus_alignments.sam

Caveats:
    - Only works for samples prepared using the linker ligation protocol.
    - Will NOT work for PCR samples (gives about 50/50 plus/minus).

Negative control (sample that should have no negative RNA):
    python strandedness.py /lab/solexa_rouskin/projects/BarCoding/Mapping_Files/ogab11_S2_L001_ogab11.sam

"""

import sys

N_FLAGS = 12  # number of flags fields in the SAM file format
FLAG_INDEX = 1  # column of SAM file in which flag is located
REFNAME_INDEX = 2  # column of SAM file in which reference name is located
POS_INDEX = 3  # column of SAM file in which the 5' position is located
TLEN_INDEX = 8  # column of SAM file in which template length is located
PRINT_EVERY = 10000  # print every time this many alignments are processed
WRITE_EVERY = 1000000  # write to file every time this many alignments are processed
write = False
pos_lo, pos_hi = None, None
"""
pos_lo, pos_hi = 28300, None
pos_lo, pos_hi = 1000, 20000
pos_lo, pos_hi = 22000, 28000
pos_lo, pos_hi = None, 250
pos_lo, pos_hi = None, 1000
"""
validate = True
try:
    sam_in = sys.argv[1]  # input SAM file to be read
    # whether to allow alignments that are discordant
    if sys.argv[2] in ['y', 'ye', 'yes']:
        allow_improper_pairs = True
    elif sys.argv[2] in ['n', 'no']:
        allow_improper_pairs = False
    else:
        raise ValueError("argument 'discordant' must be y[es] or n[o]")
    if len(sys.argv) > 3:
        try:
            # Interpret numbers as ranges of positions (optional)
            pos_lo = int(sys.argv[3])
            pos_hi = int(sys.argv[4])
        except:
            # Output SAM files to be written (optional)
            sam_outplus = sys.argv[3]
            sam_outminus = sys.argv[4]
            write = True
except:
    print(
        "usage: python strandedness.py alignment_map.sam discordant[yes/no] [new_file_for_plus_alignments.sam] [new_file_for_minus_alignments.sam]")
    exit()

# Read the input SAM file
with open(sam_in) as fi:
    # Initialize values
    read = 1  # keep track of whether the line is read 1 or read 2
    pair_count = 0  # count of total pairs (equal to count of read 1 or read 2)
    improper_pair_count = 0  # count of pairs with improper alignments
    plus_count = 0  # count of pairs aligning to plus strand
    minus_count = 0  # ditto for minus strand
    plus_strand = [None, None]  # whether read 1 and read 2 align to plus strand
    minus_strand = [None, None]  # ditto for minus strand
    improper = [None, None]  # ditto for improper alignment
    plus_lines = list()  # lines to go into the output plus strand SAM file
    minus_lines = list()  # lines to go into the output minus strand SAM file
    # Create blank output files, if writing
    if write:
        with open(sam_outplus, 'w') as fp:
            pass
        with open(sam_outminus, 'w') as fn:
            pass
    # Read the header
    line = fi.readline()
    # Header lines start with @
    while line.startswith('@'):
        if line.startswith("@SQ"):
            # Add the strandedness to the name of the reference genome
            ref_name = line.split()[1].split(':')[1]
            ref_name_plus = ref_name + "_plus"
            ref_name_minus = ref_name + "_minus"
            if write:
                plus_lines.append(line.replace(ref_name, ref_name_plus))
                minus_lines.append(line.replace(ref_name, ref_name_minus))
        else:
            if write:
                plus_lines.append(line)
                minus_lines.append(line)
        line = fi.readline()
    # Read the alignment lines
    while line != "":
        # Find the position of the alignment
        pos5p = int(line.split()[POS_INDEX])
        tlen = int(line.split()[TLEN_INDEX])
        pos3p = pos5p + abs(tlen) - 1
        if (pos_lo is not None and pos5p < pos_lo) or (pos_hi is not None and pos3p > pos_hi):
            # Skip alignments outside of the specified range
            line = fi.readline()
            continue
        # Add the strandedness to the name of the reference genome
        ref_name = line.split()[REFNAME_INDEX]
        ref_name_plus = ref_name + "_plus"
        ref_name_minus = ref_name + "_minus"
        # Find the state of each flag in the SAM file
        flags = [bit == '1' for bit in reversed(
            bin(int(line.split()[FLAG_INDEX]))[2:])]
        # Set any missing flags to False
        flags.extend([False for i in range(N_FLAGS - len(flags))])
        # Assign each flag to its own variable for convenience
        paired, proper_pair, unmap, munmap, reverse, mreverse, read1, read2, \
            secondary, qcfail, dup, supplementary = flags
        # Confirm that the read flags match the read variable
        if validate:
            assert ((read == 1) == read1 and (read == 2) == read2)
        # Say that a read is on the plus strand if it is paired and either
        # - read1 and in the reverse orientation (with read2 going forward)
        # - read2 and in the forward orientation (with read1 in reverse)
        plus_strand[read - 1] = (
                (paired and (proper_pair or allow_improper_pairs))
                and ((read1 and reverse and not mreverse) or
                     (read2 and mreverse and not reverse))
                and not (unmap or munmap or secondary or qcfail or dup
                         or supplementary)
        )
        # Say that a read is on the minus strand if it is paired and either
        # - read1 and in the forward orientation (with read2 in reverse)
        # - read2 and in the reverse orientation (with read1 going forward)
        minus_strand[read - 1] = (
                (paired and (proper_pair or allow_improper_pairs))
                and ((read1 and mreverse and not reverse) or
                     (read2 and reverse and not mreverse))
                and not (unmap or munmap or secondary or qcfail or dup
                         or supplementary)
        )
        improper[read - 1] = not proper_pair
        # Confirm that the read was not labeled as both plus and minus
        assert (not (plus_strand[read - 1] and minus_strand[read - 1]))
        if write:
            # Write the read to one of the output SAM files
            if plus_strand[read - 1]:
                plus_lines.append(line.replace(ref_name, ref_name_plus))
            elif minus_strand[read - 1]:
                minus_lines.append(line.replace(ref_name, ref_name_minus))
        if read == 2:
            # After reading the last (second) read of the alignment
            # Confirm that the paired reads agree on the following properties
            if validate:
                assert (plus_strand[0] == plus_strand[1])
                assert (minus_strand[0] == minus_strand[1])
                assert (improper[0] == improper[1])
            # Increment the counts of plus, minus, improper, and pairs
            plus_count += plus_strand[0]
            minus_count += minus_strand[0]
            improper_pair_count += improper[0]
            pair_count += 1
            # Since read is 2, switch to 1
            read = 1
        else:
            # Since read is 1, switch to 2
            read = 2
        # Report how many pairs have been processed
        if pair_count % PRINT_EVERY == 0:
            print("Total pairs:", pair_count, end="\r")
        # Write to file, if writing
        if write and pair_count % WRITE_EVERY == 0:
            with open(sam_outplus, 'a') as fp:
                fp.write(''.join(plus_lines))
            with open(sam_outminus, 'a') as fm:
                fm.write(''.join(minus_lines))
            plus_lines = list()
            minus_lines = list()
        line = fi.readline()
# Final write of unwritten lines
if write:
    with open(sam_outplus, 'a') as fp:
        fp.write(''.join(plus_lines))
    with open(sam_outminus, 'a') as fm:
        fm.write(''.join(minus_lines))

# Print the results
print("Total pairs:", pair_count)
plus_and_minus = plus_count + minus_count
print("Plus or minus pairs:", plus_and_minus)
if pair_count > 0:
    print("Plus or minus fraction:", plus_and_minus / pair_count)
# Note: the fractions are fractions of plus and minus reads, not total pairs
print("Plus pairs:", plus_count)
if plus_and_minus > 0:
    print("Plus fraction:", plus_count / plus_and_minus)
print("Minus pairs:", minus_count)
if plus_and_minus > 0:
    print("Minus fraction:", minus_count / plus_and_minus)
# Note: plus and minus can also be improper if allow_improper_pairs is True
print("Improper pairs:", improper_pair_count)
if pair_count > 0:
    print("Improper fraction:", improper_pair_count / pair_count)
