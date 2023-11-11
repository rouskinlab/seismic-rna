
.. _sam-format:

SAM, BAM, and CRAM: Alignment Maps
------------------------------------------------------------------------

Aligned reads are stored as alignment map files.

Three formats of alignment map files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These files come in three formats -- SAM, BAM, and CRAM -- documented on
the `Samtools website`_:

======================== ========= ====== ========== ======== ========= ===================
Format                   Extension Type   I/O Effort I/O Time File Size Uses in SEISMIC-RNA
======================== ========= ====== ========== ======== ========= ===================
Sequence Alignment Map   ``.sam``  text   ●          ●●       ●●●       parsing and editing
Binary Alignment Map     ``.bam``  binary ●●         ●        ●●        short-term storage
CompRessed Alignment Map ``.cram`` binary ●●●        ●●●      ●         long-term storage
======================== ========= ====== ========== ======== ========= ===================

- "I/O Effort" ranks the difficulty of reading/writing the format.
- "I/O Time" ranks the amount of time needed to read/write the format.
- "File Size" ranks the sizes of files in the format.

Throughout this documentation, we use the term "XAM format" to refer to
all three formats collectively, when the specific format -- SAM, BAM, or
CRAM -- does not matter.

Endedness of reads in XAM files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to read endedness in FASTQ files (see :ref:`fastq_endedness`),
reads in XAM files are also single- or paired-end.
SEISMIC-RNA also requires that each XAM file contain only single-end or
only paired-end reads.
Unlike with FASTQ files, paired-end XAM files must be interleaved;
SEISMIC-RNA cannot handle XAM files with only 1st or only 2nd mates.

Quality score encodings in XAM files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SAM files encode `Phred quality scores`_ in the same manner as FASTQ
files (see :ref:`phred_encodings`).
BAM and CRAM files also encode Phred scores, but in a binary format.


.. _Samtools website: https://samtools.github.io/hts-specs/
.. _Phred quality scores: https://en.wikipedia.org/wiki/Phred_quality_score
