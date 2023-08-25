
Align
========================================================================

Align performs quality control, trimming, alignment, deduplication, and
outputting of reads to a binary alignment map (BAM) file.

1.  Initial Quality Control (optional): A report of the input FASTQ file
    quality is generated using the software ``fastqc``
    (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2.  Trimming (optional): Adapters and low-quality base calls are trimmed
    from the ends of every read using the software ``cutadapt``
    (https://cutadapt.readthedocs.io/en/stable/).
3.  Trimmed Quality Control (optional): Another report of the trimmed
    FASTQ file is generated with ``fastqc``.
4.  Alignment: The trimmed FASTQ files are aligned to a FASTA file containing one or more reference sequences, yielding a sequence alignment map (SAM) file, using third-party software ``bowtie2``.
5.  Deduplication: Reads that align equally well to multiple locations in the reference sequences are removed from the SAM file using a DREEM internal function.
6. Outputting: The deduplicated SAM file is converted into a BAM file, sorted positionally, indexed, and split into one BAM file per reference sequence using the third-party software ``samtools``.


.. include:: examples.rst

.. include:: IO.rst


.. include:: cli.rst

.. include:: api.rst
