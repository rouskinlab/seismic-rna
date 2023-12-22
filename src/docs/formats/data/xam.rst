
SAM, BAM, and CRAM: Alignment Maps
--------------------------------------------------------------------------------

.. note::
    These three file formats are closely related and collectively called "XAM"
    format whenever the specific format is irrelevant.

XAM file: Content format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Comparison of SAM, BAM, and CRAM formats
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

======================== ====== ========== ======== ========= ===================
Format                   Type   I/O Effort I/O Time File Size Uses in SEISMIC-RNA
======================== ====== ========== ======== ========= ===================
Sequence Alignment Map   text   ●          ●●       ●●●       parsing and editing
Binary Alignment Map     binary ●●         ●        ●●        short-term storage
CompRessed Alignment Map binary ●●●        ●●●      ●         long-term storage
======================== ====== ========== ======== ========= ===================

- "I/O Effort" ranks the difficulty of reading/writing the format.
- "I/O Time" ranks the amount of time needed to read/write the format.
- "File Size" ranks the sizes of files in the format.

See the `Samtools website`_ for more information on SAM, BAM, and CRAM.

XAM file: Endedness
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Similar to read endedness in FASTQ files (see :ref:`fastq_endedness`), reads in
XAM files are also single- or paired-end.
SEISMIC-RNA also requires that each XAM file contain only single-end or only
paired-end reads.
Unlike with FASTQ files, paired-end XAM files must be interleaved; SEISMIC-RNA
cannot process XAM files with only 1st or only 2nd mates.

XAM file: Quality score encodings
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SAM files encode `Phred quality scores`_ in the same manner as FASTQ files (see
:ref:`phred_encodings`).
BAM and CRAM files also encode Phred scores, but in a binary format.

XAM file: Path format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

XAM file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA accepts the following extensions for XAM files:

- SAM: ``.sam``
- BAM: ``.bam``
- CRAM: ``.cram``

XAM path parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The file name is the reference.
- The file must be in a directory named ``align``.
- The directory containing ``align`` is the sample.

For example, SEISMIC-RNA would parse ``project/out/umber/align/chartreuse.cram``
to the sample ``umber`` and reference ``chartreuse``.

XAM file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

XAM as input file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Alignment maps for the Relate step must be input as XAM files.

XAM as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The Align step outputs a file in CRAM (with option ``--cram``) or BAM (with
  option ``--bam``) format for each reference to which each sample was aligned.

XAM as temporary file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The Align step writes a temporary BAM file for each input FASTQ, then splits
  that file into one BAM/CRAM file for each reference.
- The Relate step filters and converts each input BAM/CRAM file into a temporary
  SAM file, which it parses to generate the relation vectors.

.. _Samtools website: https://samtools.github.io/hts-specs/
.. _Phred quality scores: https://en.wikipedia.org/wiki/Phred_quality_score
