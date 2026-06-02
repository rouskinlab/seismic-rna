********************************************************************************
seismic align
********************************************************************************


Purpose
================================================================================

``seismic align`` maps sequencing reads to reference sequences.
It first uses fastp_ to trim low-quality bases and adapter sequences, then uses
Bowtie2_ to align reads, and finally filters out low-quality or low-coverage
alignments.
The output is one BAM file per sample/reference combination, ready for IDmut.

Requires fastp_ (≥ 0.23) and Bowtie2_ (≥ 2.5) as external dependencies.


Inputs
================================================================================

Reference FASTA
    Positional argument.

FASTQ files
    Pass reads with ``--fastqz`` (``-z``, single-end), ``--fastqy`` (``-y``,
    paired interleaved), or ``--fastqx`` (``-x``, paired separate files).
    For pre-demultiplexed FASTQ files, use ``--dmfastqz`` (``-Z``),
    ``--dmfastqy`` (``-Y``), or ``--dmfastqx`` (``-X``).


Outputs
================================================================================

All outputs go into ``{out}/{sample}/align/``.

``{ref}.bam``
    Aligned, sorted, indexed BAM file, named after the reference.

``align-report.json``
    Summary of alignment statistics.
    See :doc:`/formats/report/align`.

Unaligned reads are written to FASTQ files in the same directory
(``--bt2-un`` is on by default).


Quick example
================================================================================

Align paired-end reads in separate files to a reference::

    seismic align ref.fa -x sample1_R1.fastq.gz sample1_R2.fastq.gz


Options
================================================================================

Read trimming (fastp)
    ``--fastp/--no-fastp``
        Run fastp before alignment (default on).
    ``--fastp-3/--no-fastp-3``
        Trim low-quality bases from 3' ends of reads (default on).
    ``--fastp-5/--no-fastp-5``
        Trim low-quality bases from 5' ends of reads (default off).
    ``--fastp-adapter-trimming/--no-fastp-adapter-trimming``
        Trim adapter sequences (default on).
    ``--fastp-adapter-1 SEQ`` / ``--fastp-adapter-2 SEQ``
        Specify adapter sequences for read 1 and read 2 (default: auto-detect).
    ``--fastp-poly-g {auto|yes|no}``
        Trim poly(G) tails, common in two-color sequencing (default auto).
    ``--fastp-min-length N``
        Discard reads shorter than N bases after trimming (default 9).

Alignment (Bowtie2)
    ``--bt2-local/--bt2-end-to-end``
        Align in local mode (default) or end-to-end mode.
    ``--bt2-orient {fr|rf|ff}``
        Expected orientation of paired-end mates (default ``fr``).
    ``--bt2-X N``
        Maximum paired-end alignment length in bases (default 600).
    ``--bt2-I N``
        Minimum paired-end alignment length in bases (default 0).

Post-alignment filters
    ``--min-mapq N``
        Discard reads with mapping quality below N (default 25).
    ``--min-reads N`` (``-N``)
        Discard a BAM file if it has fewer than N reads (default 1000).

Strand separation
    ``--sep-strands/--mix-strands``
        Split each BAM into forward- and reverse-strand reads (default off).
    ``--rev-label LABEL``
        Suffix added to the reference name for reverse-strand reads
        (default ``-rev``).

Other
    ``--phred-enc N``
        Phred score encoding (default 33).
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/align_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

Very few reads align
    Check that the reference FASTA matches the organism and genome build
    of your sample.
    If fastp discards most reads, lower ``--fastp-min-length`` or disable
    adapter trimming.

BAM file skipped (too few reads)
    The BAM passed ``--min-mapq`` but fell below ``--min-reads``.
    Lower ``--min-reads`` or check for contamination.

Paired reads treated as single-end
    Verify that mate files are truly paired and named correctly when using
    ``-x``.


See also
================================================================================

- :doc:`demult` — demultiplex before aligning if needed
- :doc:`idmut` — next step: call mutations from the aligned reads
- :doc:`/formats/report/align`
- :doc:`/use/branch`, :doc:`/use/parallel`


.. _fastp: https://github.com/OpenGene/fastp
.. _Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/
