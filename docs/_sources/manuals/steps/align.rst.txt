
Align the sequencing reads
------------------------------------------------------------------------


FASTA input file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alignment requires exactly one file of reference sequences, which must
be in FASTA format.


Options for FASTQ input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SEISMIC-RNA can align single-end reads and paired-end reads where mate 1
and mate 2 reads are either interleaved in one FASTQ file or in separate
FASTQ files labeled 1 and 2.

FASTQ(s) of one sample
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Give the path to a FASTQ file with the ``-x``, ``-y``, or ``-z`` option,
depending on the contents of the FASTQ file.

FASTQ files of paired-end reads, with mates 1 and 2 in separate files::

    seismic all refs.fa -x sample_1_R1.fq.gz -x sample_1_R2.fq.gz

FASTQ file of paired-end reads, with mates 1 and 2 interleaved::

    seismic all refs.fa -y sample_1.fq.gz

FASTQ file of single-end reads::

    seismic all refs.fa -z sample_1.fq.gz

FASTQs of multiple samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Give the ``-x``, ``-y``, or ``-z`` option multiple times to process more
than one sample simultaneously.

FASTQ files of paired-end reads, with mates 1 and 2 in separate files::

    seismic all refs.fa -x sample_1_R1.fq.gz -x sample_1_R2.fq.gz -x sample_2_R1.fq.gz -x sample_2_R2.fq.gz

FASTQ files of paired-end reads, with mates 1 and 2 interleaved::

    seismic all refs.fa -y sample_1.fq.gz -y sample_2.fq.gz -y sample_3.fq.gz

FASTQ files of single-end reads::

    seismic all refs.fa -z sample_1.fq.gz -z sample_2.fq.gz -z sample_3.fq.gz

Directory of FASTQ files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If any given path is a directory instead of a FASTQ file, then it will
be searched for FASTQ files recursively (any subdirectories will also be
searched, any of their subdirectories will be searched, and so on).

.. note::
    Any file that does not have a valid FASTQ extension will be ignored,
    so if the directory contains a mixture of FASTQ and non-FASTQ files,
    then the program will still work.

Directory containing FASTQ files of paired-end reads, with mates 1 and 2
in separate files; note that if a directory contains both the mate 1 and
the mate 2 FASTQ files, then you only need to give it once::

    seismic all refs.fa -x samples/

Directory containing FASTQ files of paired-end reads, with mates 1 and 2
interleaved::

    seismic all refs.fa -y samples/

Directory containing FASTQ files of single-end reads::

    seismic all refs.fa -z samples/

Align any combination of the above
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You may any combination of the above options for specifying FASTQ files
in a single command. For example, to align

- two directories of paired-end reads with mates 1 and 2 separated
- one sample of paired-end reads with mates 1 and 2 interleaved
- two samples of single-end reads

one could type::

    seismic all refs.fa -x samples_1-10/ -x samples_11-20/ -y sample_21.fq.gz -z sample_22.fq.gz -z sample_23.fq.gz


Options for Phred score encoding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SEISMIC-RNA defaults to using Phred+33 encoding for FASTQ files, which
is standard on modern Illumina sequencers. To change the Phred score
encoding, use the option ``--phred-enc``.


Options for quality assessment with FastQC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, each FASTQ file is processed with `FastQC`_, both before and
after trimming, in order to find any potential problems. FastQC can be
disabled with the flag ``--no-fastqc``. To enable automatic extraction
of the zipped output files from FastQC, add the flag ``--qc-extract``.


Options for trimming with Cutadapt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, each FASTQ file and pair of mated FASTQ files is trimmed for
adapters and low-quality bases using `Cutadapt`_. To disable adapter
trimming, add the flag ``--no-cut``.

Adapter trimming
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, SEISMIC-RNA uses the standard, minimal adapter sequences for
Illumina sequencing runs for both read 1 and (if paired-end) read 2:

- 5': ``GCTCTTCCGATCT``
- 3': ``AGATCGGAAGAGC``

To use another adapter, type its sequence after the appropriate option:

====== ====== ==============
 Side   Read   Option
====== ====== ==============
 5'     1      ``--cut-g1``
 5'     2      ``--cut-g2``
 3'     1      ``--cut-a1``
 3'     2      ``--cut-a2``
====== ====== ==============

Quality trimming
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Base calls on either end of a read that fall below a minimum Phred score
quality are trimmed with Cutadapt. The default minimum quality is 25,
which corresponds to a confidence of 1 - 10 :sup:`-2.5` = 0.997 that the
base call is correct. To change the quality threshold, use the option
``--min-phred``.

Dark cycle trimming
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

On some Illumina sequencers (e.g. NextSeq, iSeq), the probes used to
detect G bases emit no light. Hence, these instruments will label a base
call as a G if it appears dark. If sequencing reaches the end of a read,
then there will be no more bases to sequence, so every cycle thereafter
will be dark, causing a string of Gs to be added to the 3' end of the
read. Using the option ``--cut-nextseq`` tells Cutadapt to trim off any
high-quality G bases from the 3' end of each read. This may improve the
alignment (especially in end-to-end mode) but also removes real G bases
from the 3' ends of reads (since they cannot be distinguished from any
artefactual G bases).

Additional options for Cutadapt
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

See the command line interface reference for the full list of options
that SEISMIC-RNA can use with Cutadapt. These options should suffice for
most users. If you require a more customized adapter trimming workflow,
we recommend that you perform adapter trimming and alignment outside of
SEISMIC-RNA and pass your BAM file(s) into the relate step.


Options for alignment with Bowtie 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bowtie 2 index files: automatic and pre-built
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA searches for a set of Bowtie 2 index files in the directory
where the FASTA file is. A valid Bowtie 2 index comprises six files, all
with the same name as the FASTA file, with the extensions ``.1.bt2``,
``.2.bt2``, ``.3.bt2``, ``.4.bt2``, ``.rev.1.bt2``, and ``.rev.2.bt2``.
If this set of files exists, then SEISMIC-RNA uses it as the index.
Otherwise, it calls ``bowtie2-build`` to build an index in a temporary
directory, which is deleted after alignment finishes. Indexing a small
FASTA file takes several seconds. However, indexing a very large FASTA
(e.g. a whole transcriptome) can take hours, so it is advantageous to
pre-build your index in the same directory as the FASTA file, to save
time in case you need to align more than once. You can pre-build an
index with the command ::

    bowtie2-build refs.fa refs

replacing `refs` with the path to and name of your FASTA file. See the
`Bowtie 2 Indexer manual`_ for more details.

.. note::
    If SEISMIC-RNA finds a pre-built Bowtie 2 index, then it does *not*
    verify that the index was actually built from the FASTA file of the
    same name. You must verify this yourself if using a pre-built index.
    You can assume the index is correct if you build it using the above
    command and avoid modifying or replacing the FASTA and index files.

Alignment modes: local and end-to-end
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

During alignment, Bowtie 2 can either align the entirety of each read
(end-to-end mode) or find and align only the section of the read that
yields the best alignment score (local mode). Generally, end-to-end mode
produces more spurious mutations (false positives) caused by artefacts
such as untrimmed adapters at the ends of reads; while local mode drops
more real mutations (false negatives) within several nucleotides of the
ends of reads. For mutational profiling, false positives are much more
deleterious than false negatives, so SEISMIC-RNA defaults to local mode
(``--bt2-local``). End-to-end mode (``--bt2-end-to-end``) should be used
only with amplicon-based samples and careful adapter trimming; and even
then, local mode works well enough, except for counting mutations at the
very ends of reads.

Filtering alignments
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Alignments can be filtered by `alignment score`_ and `mapping quality`_,
which are distinct properties. Alignment score measures how well a read
aligns to a given location in the reference. It is calculated from the
number of matches, substitutions, and gaps using the score parameters.
The minimum alignment scores for local and end-to-end modes can be set
using ``--bt2-score-min-loc`` and ``--bt2-score-min-e2e``, respectively.
Mapping quality measures how unique an alignment is: high quality if the
read aligns with a high score to exactly one location, low quality if it
aligns with similar scores to multiple locations in the reference. The
minimum mapping quality can be set with the option ``--min-mapq``.

Outputting unaligned reads/pairs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For some reasons, including troubleshooting low alignment rates, it can
be helpful to output all reads. The flag ``--bt2-unal`` causes all reads
(including those that did not align) to appear in the temporary SAM file
that is output directly from Bowtie 2. Because this file is located in
the temporary directory, the ``--save-temp`` flag must also be used, or
else the SAM file (and everything else in the directory) will be deleted
when the alignment step finishes.

Paired-end options
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Bowtie 2 considers paired-end reads to align "concordantly" when their
relative positions match expectations and "discordantly" otherwise. See
the `Bowtie 2 manual for details on concordant/discordant alignments`_.
By default, SEISMIC-RNA treats only concordantly aligning pairs as valid
alignments. To also treat discordant pairs as valid alignments, use the
flag ``--bt2-discordant``.

Several options control which types of alignments are concordant. First,
the expected orientation of paired mates is set using ``--bt2-orient``.
It can be ``fr`` (the 5'-most mate is forward, the 3'-most is reversed),
``rf`` (the 5'-most mate is reversed, the 3'-most is forward), or ``ff``
(both mates are forward). The default is ``fr`` (the most common type).
Second, the mates may `overlap partially or completely, or dovetail`_.
By default, overlaps (partial and complete) are considered concordant,
and dovetailing is considered discordant. The flag ``--bt2-no-contain``
treats as discordant pairs where one mate completely overlaps the other,
while ``--bt2-dovetail`` treats dovetailed pairs as concordant. Pairs
that overlap partially are always considered concordant in SEISMIC-RNA.

.. note::
    The flags ``--bt2-[no-]contain`` and ``--bt2-[no-]dovetail`` choose
    whether to treat these types of overlaps as concordant (yes) or
    discordant (no). If they are treated as discordant, then the flag
    ``--bt2-[no-]discordant`` determines whether they are considered
    valid alignments (yes) or invalid (no).

The option ``--bt2-mixed`` enables `mixed mode`_ wherein, for pairs that
fail to produce a valid paired-end alignment, Bowtie 2 attempts to align
each mate individually (as if it were a single-end read).

Troubleshooting alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alignment rate is low
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If the percent of reads aligning to the reference is less than expected,
then try the following steps (in this order):

1.  Ensure you are using Bowtie version 2.5.1 or later (version 2.5.0
    has a bug that affects alignment rate). You can check the version by
    running ``bowtie2 --version | head -n 1``.
2.  Double check that the FASTA has the correct reference sequence(s)
    and that, if the Bowtie 2 index was pre-built before the align step,
    that the correct FASTA file was used.
3.  Rerun alignment using the flags ``--bt2-unal`` and ``--save-temp``,
    which will write all the unaligned reads to the temporary SAM file
    and keep that file after alignment ends. Find unaligned reads with
    ``samtools view -f 4 temp/sample/align/align-2_align/refs.sam -o x``
    where ``sample``, ``refs``, and ``x`` are replaced with the name of
    the sample, name of the FASTA file, and name of the SAM file into
    which to write the unaligned reads, respectively. Open the SAM file,
    process several unaligned reads randomly, and use `BLAST`_ to discern
    their origins, which can help in deducing what went wrong.


.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Cutadapt: https://cutadapt.readthedocs.io/en/stable/
.. _Bowtie 2 Indexer manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
.. _alignment score: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar
.. _mapping quality: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mapping-quality-higher-more-unique
.. _Bowtie 2 manual for details on concordant/discordant alignments: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#concordant-pairs-match-pair-expectations-discordant-pairs-dont
.. _mixed mode: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mixed-mode-paired-where-possible-unpaired-otherwise
.. _overlap partially or completely, or dovetail: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
