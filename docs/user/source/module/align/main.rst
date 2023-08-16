
Alignment 
=============================

The alignment module performs quality control, trimming, alignment, deduplication, and outputting of reads to a binary alignment map (BAM) file.

1. Quality Control: A report of the FASTQ file quality is generated using the third-party software ``fastqc``.
2. Trimming: Adapters and low-quality base calls are trimmed from the ends of every read using the third-party software ``cutadapt``.
3. Trimmed Quality Control: Another report of FASTQ file quality after trimming is generated with ``fastqc``.
4. Alignment: The trimmed FASTQ files are aligned to a FASTA file containing one or more reference sequences, yielding a sequence alignment map (SAM) file, using third-party software ``bowtie2``.
5. Deduplication: Reads that align equally well to multiple locations in the reference sequences are removed from the SAM file using a DREEM internal function.
6. Outputting: The deduplicated SAM file is converted into a BAM file, sorted positionally, indexed, and split into one BAM file per reference sequence using the third-party software ``samtools``.


.. include:: examples.rst

.. include:: reference.rst

.. include:: IO.rst

.. include:: key_algos.rst

.. include:: edge_cases.rst

**Contributors**: [write who contributed to the module here.]
