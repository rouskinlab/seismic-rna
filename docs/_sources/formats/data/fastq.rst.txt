
FASTQ: Sequencing Reads
------------------------------------------------------------------------

Sequencing data files must be in `FASTQ format`_.


.. _fastq_endedness:

Endedness: single-end and paired-end reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Paired-end, two separate FASTQ files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
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

.. note::
    For every mate 1 file given, a mate 2 file with the same sample name
    must be given, and vice versa. The mates in file 2 must be in the
    same order as those in file 1.

Paired-end, one interleaved FASTQ file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

Single-end FASTQ file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

FASTQ file naming
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FASTQ file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The following file extensions are accepted for FASTQ files:

- Compressed (`gzip`_): ``.fastq.gz``, ``.fq.gz``
- Uncompressed: ``.fastq``, ``.fq``

.. note::
    SEISMIC-RNA accepts FASTQ files that are compressed with `gzip`_.
    It is recommended to always use compressed FASTQ files unless they
    must be read by a person because FASTQ files are typically very
    large, on the order of 100 Mb to 10 Gb. The file extension will be
    preserved through the pipeline: if an input FASTQ file has the
    extension ``.fq.gz``, then the trimmed FASTQ file (if any) will also
    have that extension and be compressed.

FASTQ mate 1 and 2 labels
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For paired-end reads in which mate 1s and mate 2s are in separate files,
the file names must have one of the following labels right before the
file extension:

- Mate 1: ``_R1``, ``_mate1``, ``_1_sequence``, ``_R1_001``, ``_mate1_001``, ``_1_sequence_001``
- Mate 2: ``_R2``, ``_mate2``, ``_2_sequence``, ``_R2_001``, ``_mate2_001``, ``_2_sequence_001``

If you would like future versions to support additional file extensions,
then please request so by creating an issue (see ).

FASTQ name parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For FASTQ files given via ``--fastqp`` (``-x``), ``--fastqi`` (``-y``),
or ``--fastqs`` (``-z``), the sample name is determined by parsing the
FASTQ file name. For demultiplexed FASTQ files given via ``--dmfastqp``
(``-X``), ``--dmfastqi`` (``-Y``), or ``--dmfastqs`` (``-Z``), the
reference name is determined by parsing the FASTQ file name, and the
sample name comes from the directory in which the FASTQ file is located.

When parsing the name of the sample/reference from the FASTQ file name,
the name of the file up to but not including the file extension is used
for single-end and interleaved paired-end files. For paired-end reads
in separate files, the mate 1 and 2 labels are removed first.

FASTQ symbols
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FASTQ DNA alphabet
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Read sequences may contain any uppercase characters, but all characters
besides A, C, G, and T are treated as any nucleotide (i.e. N).

.. _phred_encodings:

Phred quality score encodings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

`Phred quality scores`_ represent the confidence that a base in a FASTQ
file was called correctly during sequencing.
The probability *p* that a base was called incorrectly is 10 raised to
the power of the quotient of the Phred score *s* and -10:

*p* = 10 :sup:`-s/10`

For example, if a base call has a Phred score of 30, the probability
that the base call is incorrect is 10 :sup:`-30/10` = 0.001.

In FASTQ files, each phred quality score (a non-negative integer) is
encoded as one character of text by adding another integer *N* to the
Phred score (`Phred+N`_) and then converting the number to the character
with the corresponding `ASCII code`_.
For example, if *N* is 33, then the Phred score 25 would be encoded by
adding 33 to 25 (obtaining 58), then writing the character whose ASCII
code is 58 (which is ``:``).


.. _FASTQ format: https://en.wikipedia.org/wiki/FASTQ_format
.. _guide from Illumina: https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html
.. _gzip: https://www.gnu.org/software/gzip/
.. _Phred quality scores: https://en.wikipedia.org/wiki/Phred_quality_score
.. _Phred+N: https://en.wikipedia.org/wiki/FASTQ_format#Encoding
.. _ASCII code: https://en.wikipedia.org/wiki/ASCII
