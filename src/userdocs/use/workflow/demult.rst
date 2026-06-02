********************************************************************************
seismic demult
********************************************************************************


Purpose
================================================================================

``seismic demult`` splits one or more multiplexed FASTQ files into per-sample
FASTQ files by matching barcode sequences embedded in the reads.
Use it when multiple samples were pooled before sequencing and share a single
FASTQ file, distinguished only by their barcodes.

Requires seqkit_ (≥ 2.10.1).


Inputs
================================================================================

Multiplexed FASTQ files
    Pass reads with ``--fastqz`` (``-z``, single-end), ``--fastqy`` (``-y``,
    paired interleaved), or ``--fastqx`` (``-x``, paired separate files).
    For FASTQ files that are already partially demultiplexed and need re-sorting,
    use ``--dmfastqz`` (``-Z``), ``--dmfastqy`` (``-Y``), or
    ``--dmfastqx`` (``-X``).

Reference FASTA
    Positional argument.
    Required so SEISMIC-RNA knows which reference name corresponds to each
    barcode.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/demult/``.
Reads that do not match any barcode go into an ``unmatched`` subdirectory.


Quick example
================================================================================

Demultiplex a single-end FASTQ with two barcodes defined on the command line::

    seismic demult ref.fa -z reads.fastq.gz \
        --barcode sample1 ACGTACGT 5 \
        --barcode sample2 TTGGCCAA 5


Options
================================================================================

Define barcodes
    ``--barcode NAME SEQ POS``
        Specify one barcode: its name (must match a reference in the FASTA),
        DNA sequence, and 1-based position in the read.
        Repeat for each sample.
    ``--refs-meta FILE``
        Alternatively, read barcode definitions from a CSV file.
        See :doc:`/formats/meta/refs`.
    ``--barcode-start N`` / ``--barcode-end N``
        When all barcodes occupy the same region, set the 0-based start and end
        positions (half-open) of that region (default 0 for both).
    ``--read-pos N``
        Expected 1-based position of the barcode in the read; defaults to
        ``--barcode-start`` when not set.

Matching tolerance
    ``--mismatch-tolerance N``
        Allow up to N mismatches between a read and the expected barcode
        (default 0).
        Use caution above 2 — computation grows factorially with mismatches.
    ``--index-tolerance N``
        Allow the barcode to appear up to N bases away from the expected position
        (default 0).
    ``--allow-n/--no-allow-n``
        Count N bases as valid mismatches (default off).
        Requires ``--mismatch-tolerance ≥ 1``; increases memory use.

Other
    ``--phred-enc N``
        Phred score encoding of the FASTQ files (default 33).
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/demult_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`align` — next step: map demultiplexed reads to the reference
- :doc:`/formats/meta/refs` — CSV format for ``--refs-meta``
- :doc:`/use/branch`, :doc:`/use/parallel`


.. _seqkit: https://bioinf.shenwei.me/seqkit/
