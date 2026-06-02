********************************************************************************
seismic splitbam
********************************************************************************


Purpose
================================================================================

``seismic splitbam`` splits a multi-reference BAM file into one BAM file per
reference sequence.
Use it when you aligned reads to a multi-reference FASTA using an external
aligner and need separate per-reference BAMs before running ``seismic idmut``.
(``seismic align`` produces per-reference BAMs automatically; this step is
only needed for externally produced BAMs.)

Requires Bowtie2_ (≥ 2.5) and SAMtools_ as external dependencies.


Inputs
================================================================================

Reference FASTA file
    First positional argument.

BAM / SAM / CRAM files
    One or more alignment files containing reads mapped to multiple references.
    See :doc:`/use/inputs`.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/align/``.

``{ref}.bam``
    One BAM file per reference, sorted and indexed.

``split-report.json``
    Summary of settings and results.


Quick example
================================================================================

Split a multi-reference BAM::

    seismic splitbam ref.fa multi.bam


Options
================================================================================

Alignment filters
    ``--min-mapq N``
        Discard reads with mapping quality below N (default 25).
    ``--min-reads N`` (``-N``)
        Discard BAM files with fewer than N reads (default 1000).
    ``--bt2-orient {fr|rf|ff}``
        Expected paired-end mate orientation (default ``fr``).
    ``--bt2-X N``
        Maximum paired-end alignment length (default 600).

Strand separation
    ``--sep-strands/--mix-strands``
        Split each BAM into forward- and reverse-strand reads (default off).
    ``--rev-label LABEL``
        Suffix for reverse-strand reference names (default ``-rev``).

Branches
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/align_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`/use/workflow/align` — produces per-reference BAMs automatically
- :doc:`/use/workflow/idmut` — consumes per-reference BAMs
- :doc:`/formats/data/xam`
- :doc:`/use/branch`, :doc:`/use/parallel`


.. _Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/
.. _SAMtools: https://www.htslib.org/
