
Relate each read to every reference position
------------------------------------------------------------------------

Input files for relate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Relate requires exactly one FASTA file containing one or more reference
sequences and any number of alignment map files in BAM format.

.. note::
    The references in the FASTA file must match those to which the reads
    were aligned to produce the BAM file(s); the names and the sequences
    must be identical. If the names differ, then the BAM files using the
    old names will be ignored; while if the sequences differ, then reads
    can yield erroneous relation vectors or fail outright.

One BAM file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A single BAM file (``sample_1/align/ref_1.bam``) can be run as follows::

    seismic relate refs.fa sample_1/align/ref_1.bam

Multiple BAM files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Multiple BAM files can be run in parallel by giving multiple paths::

    seismic relate refs.fa sample_1/align/ref_1.bam

and/or by using `glob patterns`_::

    seismic relate refs.fa sample_*/align/ref_1.bam sample_*/align/ref_2.bam

BAM file content conventions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Generally, a BAM file can contain reads that have aligned to any number
of reference, as well as unaligned reads. However, SEISMIC-RNA requires
that each BAM file contain reads aligned to exactly one reference. This
restriction enables the relate step to process BAM files in parallel,
which increases the speed. If the BAM files were created using ``seismic
align``, then they will already follow this convention.

BAM file path conventions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The name of the BAM file (minus the extension ``.bam``) must be the name
of the reference to which it was aligned. It must be inside a directory
named ``align``, which must be inside a directory named after the sample
from which the reads came. If the BAM files were created using ``seismic
align``, then they will already follow this convention.
