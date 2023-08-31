
FASTQ: Sequencing Reads
------------------------------------------------------------------------

Sequencing data files must be in `FASTQ format`_. Single- and paired-end
reads are supported, and paired-end reads can be provided either as two
files where mate 1 and mate 2 reads are separated, or as one interleaved
file where the mate 1 and mate 2 reads alternate. However, single- and
paired-end reads cannot be mixed in one FASTQ file, nor can separate and
interleaved reads.

.. _FASTQ format: https://en.wikipedia.org/wiki/FASTQ_format

Endedness of sequencing reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

If you would like future versions of DREEM to support additional file
extensions, please create a new issue on GitHub. [REF]

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

FASTQ quality score encodings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The `Phred quality scores`_ are encoded by adding an integer *N* to the
Phred score (`Phred+N`_). Most modern Illumina instruments output FASTQ
files with Phred+33 encoding (which is the default in SEISMIC-RNA), but
Phred+64 is also common. The quality score encoding can be changed (in
this example, to Phred+64) with the option ``--phred-enc``::

    seismic align --phred-enc 64


.. _gzip: https://www.gnu.org/software/gzip/
.. _Phred quality scores: https://en.wikipedia.org/wiki/Phred_quality_score
.. _Phred+N: https://en.wikipedia.org/wiki/FASTQ_format#Encoding
