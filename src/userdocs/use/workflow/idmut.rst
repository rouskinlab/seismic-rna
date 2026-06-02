********************************************************************************
seismic idmut
********************************************************************************

.. note::

    ``idmut`` was named ``relate`` in versions before v0.25.0.
    See the project ``CHANGELOG.md`` for the full rename list.


Purpose
================================================================================

``seismic idmut`` ("identify mutations") reads each aligned read and records,
position by position, whether each base was a match, a substitution (and which
one), a deletion, an insertion, or ambiguous.
Its output is the foundation for every downstream analysis in SEISMIC-RNA.


Inputs
================================================================================

Reference FASTA file
    First positional argument.
    Must be the same FASTA used during alignment.
    See :doc:`/formats/data/fasta`.

Alignment map files (SAM / BAM / CRAM)
    One or more alignment files as additional positional arguments.
    Each file must be named after the single reference it was aligned to
    (minus the extension).
    See :doc:`/use/inputs`.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/idmut/{ref}/``.

``idmut-batch-{num}.brickle``
    One batch of per-read, per-position mutation calls.
    See :doc:`/formats/data/brickle`.

``idmut-report.json``
    Summary of settings and results.
    See :doc:`/formats/report/idmut`.


Quick example
================================================================================

Identify mutations in two alignment maps from ``sample-1``::

    seismic idmut refs.fa out/sample-1/align/ref-1.bam out/sample-1/align/ref-2.bam


Options
================================================================================

Quality filtering
    ``--min-phred N``
        Flag base calls with Phred score below N as ambiguous (default 25).
    ``--min-mapq N``
        Discard reads with mapping quality below N (default 25).
    ``--min-reads N`` (``-N``)
        Skip alignment files with fewer than N reads (default 1000).

Indel handling
    ``--ambindel/--no-ambindel``
        Report all possible placements of each indel in repetitive regions
        (default on).
        Disable with ``--no-ambindel`` for a speed-up with minimal accuracy
        loss on low-indel data.
    ``--insert3/--insert5``
        Anchor each insertion on the base to its 3' side (``--insert3``, the
        default) or its 5' side (``--insert5``).

Read trimming
    ``--clip-end5 N`` / ``--clip-end3 N``
        Ignore N positions from the 5' / 3' end of every read (default 4 each).
    ``--overhangs/--no-overhangs``
        Include bases in mate overhangs when computing mutations (default on).

Strand handling
    ``--sep-strands/--mix-strands``
        Split output into forward- and reverse-strand reads (default off).
    ``--rev-label LABEL``
        Suffix added to reverse-strand reference names (default ``-rev``).

Branches
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/idmut_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--batch-size N``
        Maximum reads per batch (default 65536).
        Smaller batches use less memory but create more output files.
        See :doc:`/use/parallel`.

    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

IDmut crashes or uses too much memory
    Reduce ``--batch-size`` and/or ``--num-cpus``.

No output for a reference
    Either every alignment file for that reference was filtered by
    ``--min-reads`` or ``--min-mapq``, or the file names did not match
    the reference names.

Indels spread across repetitive regions
    Expected behavior with ``--ambindel`` on.
    See :doc:`/algos/index`.


See also
================================================================================

- :doc:`align` — produces the alignment files this step consumes
- :doc:`filter` — next step: filter reads and positions
- :doc:`/formats/data/xam`, :doc:`/formats/data/brickle`,
  :doc:`/formats/report/idmut`
- :doc:`/use/parallel`, :doc:`/use/inputs`, :doc:`/use/branch`
