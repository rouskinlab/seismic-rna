
Demultiplex: Split multiplexed FASTQ files by their barcodes
--------------------------------------------------------------------------------

The Demultiplex step splits one or more multiplexed FASTQ files into per-sample
FASTQ files by matching barcode sequences in the reads.
It requires seqkit_ (≥ 2.10.1) as an external dependency.

Demultiplex: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Demultiplex input file: FASTQ files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Provide multiplexed FASTQ files using the same flags as the Align step:

- ``-z`` / ``--fastqz``: single-end FASTQ files
- ``-y`` / ``--fastqy``: paired-end FASTQ files with mates interleaved in one file
- ``-x`` / ``--fastqx``: paired-end FASTQ files with mates in separate files

You may also supply pre-demultiplexed FASTQ files that need to be re-sorted:

- ``-Z`` / ``--dmfastqz``, ``-Y`` / ``--dmfastqy``, ``-X`` / ``--dmfastqx``

Demultiplex input file: Reference FASTA
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Provide a reference FASTA file so that SEISMIC-RNA knows which reference
sequence corresponds to each barcode.

Demultiplex: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Demultiplex setting: Define barcodes
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can define barcodes in two ways:

**1. Per-barcode on the command line** using ``--barcode``:

Each ``--barcode`` takes three values: the barcode name (which must match a
reference name in your FASTA), the barcode DNA sequence, and the 1-indexed
position of the barcode in the read.
You can repeat ``--barcode`` for every sample in the pool::

    seismic demult ref.fa reads.fastq.gz \
        --barcode sample1 ACGTACGT 5 \
        --barcode sample2 TTGGCCAA 5

**2. Via a metadata CSV file** using ``--refs-meta``:

The CSV file should have columns for the reference name, barcode sequence, and
barcode position.
See :doc:`../../formats/meta/refs` for the format specification.

Demultiplex setting: Barcode coordinates
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When all barcodes are at the same position in the read, you can use
``--barcode-start`` and ``--barcode-end`` to set the range of positions
(0-indexed, half-open) that define the barcode in the read, rather than
specifying the position for each barcode individually.

``--read-pos`` sets the expected 1-indexed position of the barcode within the
read for all barcodes; it defaults to ``--barcode-start`` when not specified.

Demultiplex setting: Mismatch tolerance
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, barcodes must match exactly.
Use ``--mismatch-tolerance`` to allow up to this many mismatches between a read
sequence and the expected barcode.

.. warning::

    Increasing ``--mismatch-tolerance`` raises the risk of mis-assigning reads
    and increases computation time at a factorial rate; use caution above 2.

Use ``--index-tolerance`` to allow the barcode to be found up to this many
bases away from the expected position in the read.

Use ``--allow-n`` to count an N base in the read as a valid mismatch (rather
than an automatic mismatch) when ``--mismatch-tolerance`` ≥ 1.
This increases memory consumption but can recover reads with low-quality
positions in the barcode region.

Demultiplex: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output FASTQ files go into the directory
``{out}/{sample}/demult/``,
where ``{out}`` is the output directory and ``{sample}`` is the sample name
that matched the barcode.
Reads that did not match any barcode are written to an ``unmatched`` directory.

.. _seqkit: https://bioinf.shenwei.me/seqkit/
