
FASTQ: Sequencing Reads
------------------------------------------------------------------------

FASTQ file: Content format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a FASTQ file, each sequencing read record comprises four lines:

1. A header that begins with ``@`` and contains the name of the read.
2. The sequence of the read.
3. A second header that begins with ``+`` and may repeat the read name.
4. The quality string of the same length as the read, indicating the
   quality of each base in the read using :ref:`phred_encodings`.

See `FASTQ format`_ for more information.

FASTQ character sets
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

FASTQ DNA alphabet
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Read sequences may contain only ``A``, ``C``, ``G``, ``T``, and ``N``.

.. _phred_encodings:

Phred quality score encodings
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

`Phred quality scores`_ represent the confidence that a base in a FASTQ
file was called correctly during sequencing.
The probability *p* that a base was called incorrectly is 10 raised to
the power of the quotient of the Phred score *s* and -10:

*p* = 10\ :sup:`-s/10`

For example, if a base call has a Phred score of 30, the probability
that the base call is incorrect is 10\ :sup:`-30/10` = 0.001.

In FASTQ files, each phred quality score (a non-negative integer) is
encoded as one character of text by adding another integer *N* to the
Phred score (`Phred+N`_) and then converting the number to the character
with the corresponding `ASCII code`_.
For example, if *N* is 33, then the Phred score 25 would be encoded by
adding 33 to 25 (obtaining 58), then writing the character whose ASCII
code is 58 (which is ``:``).

.. _fastq_endedness:

Endedness: single-end and paired-end reads
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Illumina sequencers can be run in two main modes: single-end mode, in
which each molecule is read starting from one end; and paired-end mode,
in which each molecule is read from both ends.
For more details on endedness, see this `guide from Illumina`_.

Each paired-end read comprises a pair of so-called "mates": mate 1 and
mate 2.
Each mate is similar to a single-end read, the difference being that it
is paired with and shares its name with its mate.
There are two options for paired-end reads:

- A single "interleaved" file in which mates 1 and 2 alternate, each
  mate 2 coming directly after the mate 1 to which it is paired.
- Two separate files containing all mate 1 reads and all mate 2 reads,
  respectively, in the same order in both files.

In SEISMIC-RNA (as well as many other pieces of software), each FASTQ
file must contain only one type of reads:

- single-end reads
- paired-end reads, with interleaved 1st and 2nd mates
- paired-end reads, with 1st mates only
- paired-end reads, with 2nd mates only

Thus, we refer to an entire FASTQ file as "single-end" if it contains
only single-end reads, "interleaved paired-end" if it contains only
interleaved 1st and 2nd mates, and "separate paired-end" if it contains
either 1st or 2nd mates of paired-end reads.

The following diagrams illustrate single-end, separate paired-end, and
interleaved paired-end FASTQ files.

FASTQ file with single-end reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

1 file, 4 lines per read

``name.fq`` ::

    @Read_1_ID
    GCATGCTAGCCA
    +
    FFFFFFFFF:F:
    @Read_2_ID
    ATCGTCATGTGT
    +
    FFFFFFF:FFFF

FASTQ file with interleaved paired-end reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

1 file, 8 lines per paired read (4 for mate 1, then 4 for mate 2)

``name.fq`` ::

    @Read_1_ID/1
    GCATGCTAGCCA
    +
    FFFFFFFFF:F:
    @Read_1_ID/2
    TACGTCGTCGTC
    +
    FFFFF:FF:F::
    @Read_2_ID/1
    ATCGTCATGTGT
    +
    FFFFFFF:FFFF
    @Read_2_ID/2
    CACGAGCGATAG
    +
    FFFF:FF:::F:

Pair of FASTQ files with 1st and 2nd mates in separate files
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
2 files, 4 lines per mate in each file

``name_R1.fq`` ::

    @Read_1_ID/1
    GCATGCTAGCCA
    +
    FFFFFFFFF:F:
    @Read_2_ID/1
    ATCGTCATGTGT
    +
    FFFFFFF:FFFF

``name_R2.fq`` ::

    @Read_1_ID/2
    TACGTCGTCGTC
    +
    FFFFF:FF:F::
    @Read_2_ID/2
    CACGAGCGATAG
    +
    FFFF:FF:::F:

FASTQ file: Path format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FASTQ file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA accepts the following extensions for FASTQ files:

- Compressed (with `gzip`_):

  - ``.fq.gz`` (default)
  - ``.fastq.gz``

- Uncompressed:

  - ``.fq`` (default)
  - ``.fastq``

.. note::
    SEISMIC-RNA accepts FASTQ files that are compressed with `gzip`_.
    It is recommended to always use compressed FASTQ files because FASTQ
    files are typically very large without compression, on the order of
    100 Mb to 10 Gb.
    The file extension will be preserved through the workflow, i.e. if
    an input FASTQ file has the extension ``.fq.gz``, then the trimmed
    FASTQ file (if any) will also have that extension and be compressed.

FASTQ mate 1 and 2 labels
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For paired-end reads whose 1st and 2nd mates are in separate files, the
file names must have one of the following labels before the extension:

- Mate 1: ``_R1``, ``_mate1``, ``_1_sequence``, ``_R1_001``, ``_mate1_001``, ``_1_sequence_001``
- Mate 2: ``_R2``, ``_mate2``, ``_2_sequence``, ``_R2_001``, ``_mate2_001``, ``_2_sequence_001``

For example, a sample named ``sample-26`` consisting of paired-end reads
could have the FASTQ files ``sample-26_R1.fq`` and ``sample-26_R2.fq``.

If you would like future versions to support additional file extensions,
then please request so by creating an issue (see :doc:`../../issues`).

FASTQ path parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For FASTQ files from whole samples (``-x``, ``-y``, ``-z``), the sample
name is taken from the file name (dropping the mate number, if any).
For example, the single-end FASTQ ``-z project/sienna.fq``
would be parsed to have sample name ``sienna``.
And the separate paired-end FASTQ ``-x project/lavender_R1.fq``
would be parsed to have sample name ``lavender``.

For demultiplexed FASTQ files (``-X``, ``-Y``, ``-Z``), the reference
name is taken from the file name (dropping the mate number, if any),
and the sample name is taken from the directory of the FASTQ file.
For example, the single-end FASTQ ``-Z project/azure/ochre.fq``
would be parsed to have sample ``azure`` and reference ``ochre``.
And the separate paired-end FASTQ ``-X project/lilac/teal_R2.fq``
would be parsed to have sample ``lilac`` and reference ``teal``.

FASTQ file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FASTQ as input file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sequencing reads for these commands must be input as FASTQ files:

- ``all``
- ``align``

FASTQ as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``align`` command outputs a file in FASTQ format containing the
  unaligned reads from each input FASTQ (with option ``--bt2-un``).

FASTQ as temporary file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``align`` command writes a temporary FASTQ file for each input
  FASTQ that it trims with cutadapt (with option ``--cut``).

.. _FASTQ format: https://en.wikipedia.org/wiki/FASTQ_format
.. _guide from Illumina: https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html
.. _gzip: https://www.gnu.org/software/gzip/
.. _Phred quality scores: https://en.wikipedia.org/wiki/Phred_quality_score
.. _Phred+N: https://en.wikipedia.org/wiki/FASTQ_format#Encoding
.. _ASCII code: https://en.wikipedia.org/wiki/ASCII
