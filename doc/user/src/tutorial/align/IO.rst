
I/O files
++++++++++++++++++++++++

Input files
-----------

The alignment module requires two input files:

- A list of sequencing reads in `FASTQ format <https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq-files>`_
- A list of one or more reference sequences in `FASTA format <https://www.ncbi.nlm.nih.gov/genbank/fastaformat/>`_

Endedness: single vs. paired
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The alignment module accepts both single-end and paired-end reads.
Each FASTQ file must contain one type (not a mix of both), although if multiple FASTQ files are given, they may be of different types from each other.
The accepted types of FASTQ files are as follows:

- **Single-end**: 1 file, 4 lines per read. Example showing 2 reads:

``name.fq`` ::

    @Read_1_ID
    GCATGCTAGCCA
    +
    FFFFFFFFF:F:
    @Read_2_ID
    ATCGTCATGTGT
    +
    FFFFFFF:FFFF

- **Paired-end, separate** files for mate 1 and mate 2 reads: 2 files, 4 lines per mate in each file; mates must be in the same order in both files. Example showing 2 reads:

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

Note that for every mate 1 file given, a mate 2 file with the same sample name must be given, and vice versa.

- **Paired-end, interleaved** mate 1 and mate 2 reads: 1 file, 8 lines per paired read (4 for mate 1, then 4 for mate 2). Example showing 2 reads:

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


File extensions
~~~~~~~~~~~~~~~

The following file extensions are accepted for each file format:

- FASTA: ``.fasta``, ``.fna``, ``.fa``
- FASTQ: ``.fastq``, ``.fq``, ``.fastq.gz``, ``.fq.gz``

Note that FASTQ files may be compressed with ``gzip`` (and end with ``.gz``), while FASTA files may not.
As FASTQ files from deep sequencing may occupy many gigabytes when uncompressed, it is recommended to keep them compressed (``.gz``) unless they must be demultiplexed.
The file extension will be preserved throughout the pipeline; that is, if an input FASTQ file is compressed and has the extension ``.fq.gz``, then the temporary trimmed FASTQ file will also be compressed and end in ``.fq.gz``.

For paired-end reads in which mate 1s and mate 2s are in separate files, the file names must have one of the following suffixes immediately before the file extension:

- FASTQ 1: ``_R1``, ``_mate1``, ``_1_sequence``, ``_R1_001``, ``_mate1_001``, ``_1_sequence_001``
- FASTQ 2: ``_R2``, ``_mate2``, ``_2_sequence``, ``_R2_001``, ``_mate2_001``, ``_2_sequence_001``

If you would like future versions of DREEM to support additional file extensions, please `create a new issue on GitHub <https://github.com/rouskinlab/dreem/issues>`_.

File nomenclature
~~~~~~~~~~~~~~~~~

The paths of the input files are parsed to determine the names of the sample and reference.

- For a FASTQ file containing reads from an entire sample (arguments ``fastqs``, ``fastqi``, or ``fastq1`` and ``fastq2``), the sample name is the name of the file, minus its extension (e.g. ``.fq``) and suffix indicating mate 1 or 2 (e.g. ``_R1``). For example, if the path is ``/home/dms/foo_R2.fq.gz``, then the sample is ``foo``. The name of the reference sequence comes from the label of the sequence inside the FASTA file (the name of the FASTA file itself is used only in temporary files).

- For a FASTQ file containing reads from only one reference or construct (i.e. output from the demultiplexing module), the sample name is the directory in which the FASTQ file is located, and the reference is the name of the file, minus the suffix and extension. For example, if the path is ``/home/dms/bar/baz_mate1.fq.gz``, then the sample is ``bar`` and the reference is ``baz``. (A sequence named ``baz`` must then be present in the given FASTA file.)


Sequence alphabets
~~~~~~~~~~~~~~~~~~

Reference sequences must contain only the uppercase characters ``A``, ``C``, ``G``, and ``T``.
Read sequences may contain any uppercase characters, but all characters besides ``A``, ``C``, ``G``, and ``T`` (including `degenerate bases defined by the IUPAC <https://en.wikipedia.org/wiki/Nucleic_acid_notation>`_) are treated as any nucleotide (i.e. ``N``).

Quality score encodings
~~~~~~~~~~~~~~~~~~~~~~~

The `Phred quality scores <https://en.wikipedia.org/wiki/Phred_quality_score>`_ in FASTQ files are encoded by adding an integer *N* to the Phred score `(Phred+N) <https://en.wikipedia.org/wiki/FASTQ_format#Encoding>`_.
Most modern Illumina instruments output FASTQ files with Phred+33 encoding (which is the default in DREEM), but Phred+64 is also common.
The quality score encoding can be set to a non-default value (in this example, Phred+64) as follows:

- CLI: ``--phred-enc 64``
- API: ``phred_enc=64``


Output files
------------

All output files and directories are written in the user-specified directory ``output_dir``.

Alignment maps
~~~~~~~~~~~~~~

For each sample named ``sample_name``, given as a FASTQ file:
    For each reference sequence named ``reference_name`` in the input FASTA file:
        - A set of all reads that aligned to the reference in `binary alignment map (BAM) format <https://samtools.github.io/hts-specs/>`_: ``{output_dir}/alignment/{sample_name}/{reference_name}.bam``

FASTQC reports
~~~~~~~~~~~~~~

For each input FASTQ file named ``file_name.x`` coming from sample ``sample_name``:
    - A FASTQC report of the input file. Zipped by default; can extract automatically using ``--qc-extract`` (CLI) or ``qc_extract=True`` (API): ``{output_dir}/alignment/{sample_name}/qc-inp/{file_name}_fastqc.zip``
    - A FASTQC report of the file after trimming with Cutadapt: ``{output_dir}/alignment/{sample_name}/qc-cut/{file_name}_fastqc.zip``
