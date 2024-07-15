
FASTA: Reference sequences
------------------------------------------------------------------------

FASTA file: Content format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a FASTA file, each sequence record comprises two parts:

- A header that contains the name of the sequence.
- A body that contains the sequence.

See `FASTA format`_ for more information.

FASTA header lines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Every header line must start with the character ``>``.
The name of the sequence must follow this character, on the same line.
Optionally, metadata may follow the name of the sequence after a break
by non-alphanumeric characters such as whitespace or ``|``.

SEISMIC-RNA requires that header lines contain no metadata -- i.e. that
all characters after the initial ``>`` are part of the sequence name.
This restriction exists because the name of each reference sequence is
incorporated into file paths, and SEISMIC-RNA restricts the characters
allowed in file paths to avoid any potential problems caused by special
characters and whitespace in paths.
If SEISMIC-RNA were to simply ignore characters after the first non-path
character in the header of a FASTA, then the names would not necessarily
match those produced by other tools such as Bowtie 2 that read the FASTA
files directly; and these inconsistencies in names could cause errors.
Thus, to ensure consistent names, SEISMIC-RNA will raise errors if a
FASTA file has any illegal characters in its header lines.

FASTA body lines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The remaining lines in each record encode the sequence.
SEISMIC-RNA can parse sequences that obey the following rules:

- Alphabet: ``A``, ``C``, ``G``, and ``N`` are valid characters for DNA
  and RNA; ``T`` and ``U`` are also valid for DNA and RNA, respectively.
  Lowercase equivalents are also valid but will be cast to uppercase.
  All other characters (including whitespace) are illegal.
- Sequence lengths: Arbitrary lengths are supported, from zero to the
  maximum number of nucleotides that will fit in your system's memory.
- Line lengths: Arbitrary lengths are supported, up to the line length
  limit imposed by your system.
- Blank lines: Blank lines (i.e. containing only a newline character)
  are simply ignored, but lines containing other whitespace characters
  are illegal.

FASTA file: Path format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FASTA file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA accepts the following extensions for FASTA files:

- ``.fa`` (default)
- ``.fna``
- ``.fasta``

FASTA path parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The name of an input FASTA file of all reference sequences is used for
the following purposes:

- Determining if a Bowtie 2 index exists for the FASTA file.
- Building a Bowtie 2 index for the FASTA file.
- Linking a CRAM file to its reference sequence.

FASTA file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FASTA as input file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Reference sequences for these commands must be input as FASTA files:

- ``seismic wf``
- ``seismic align``
- ``seismic relate``
- ``seismic fold``

FASTA as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``align`` command outputs a file in FASTA format alongside each
  file in CRAM format (with option ``--cram``).

FASTA as temporary file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``align`` command writes a temporary FASTA file with a single
  reference sequence for each demultiplexed FASTQ file, which is used to
  build an index for Bowtie 2.
- The ``fold`` command writes a temporary FASTA file, which is used by
  the program ``Fold``.

.. _FASTA format: https://en.wikipedia.org/wiki/FASTA_format
