
Key algorithm(s)
++++++++++++++++

Trimming reads
--------------
DREEM uses `Cutadapt <https://cutadapt.readthedocs.io/en/stable>`_ to trim adapters and low-quality bases from the ends of reads.

Aligning reads
--------------
DREEM uses `Bowtie2 <https://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ to align reads to sets of reference sequences.

Deduplicating alignments
------------------------
DREEM uses a custom algorithm to remove reads and pairs that align equally well to two or more locations in the set of references.
For each alignment, DREEM finds the alignment score from the field ``AS:i:`` and the score of the best other alignment of the same read (if any) from the field ``XS:i:``.
The alignment is the best for the read if either no other alignment exists or the alignment score is strictly greater than the score of the best other alignment for the read.
For paired-end reads, if both mates are present, then one paired alignment must also be the best-scoring alignment for both reads individually in order to count as the best alignment.

Manipulating SAM and BAM files
------------------------------
DREEM uses `Samtools <https://www.htslib.org/>`_ to perform standard format conversion, filtering, and indexing operations on SAM and BAM files
(not including the custom deduplicating algorithm mentioned above, and the batching and vectoring algorithms described in the vector module).
