
Align: Trim FASTQ files and align them to reference sequences
--------------------------------------------------------------------------------

Align: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Align input file: Reference sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You need one file of reference sequences in FASTA format (for details on this
format, see :doc:`../../formats/data/fasta`).
If your file has characters or formatting incompatible with SEISMIC-RNA, then
you can fix it using the :doc:`../cleanfa` tool.

Align input file: Read sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Your read sequences must be in FASTQ format (see :doc:`../../formats/data/fastq`
for details on this format).

Reads can be single-end or paired-end
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can align FASTQ files of single- and paired-end reads with SEISMIC-RNA.
(For definitions of single- and paired-end reads, see :ref:`fastq_endedness`).
SEISMIC-RNA requires that single- and paired-end reads not be mixed within one
file, but it can accept different types of reads in separate FASTQ files.

Reads can come from whole or demultiplexed samples
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can align FASTQ files that come from a whole sample (possibly containing
multiple RNA sequences) or that have been demultiplexed before alignment so
that they contain reads from only one RNA sequence.
For more information on demultiplexing and how to perform it if needed, see
:doc:`./demult`.

How to specify the endedness and source of reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Specify the endedness (single-end, paired-end interleaved, or paired-end in
separate files) and source of the reads (whole or demultiplexed sample) in a
FASTQ file by giving the file path after the appropriate option:

========================== ===================== =======================
Endedness                  Whole Sample          Demultiplexed
========================== ===================== =======================
paired-end, separate files ``--fastqx`` (``-x``) ``--dmfastqx`` (``-X``)
paired-end, interleaved    ``--fastqy`` (``-y``) ``--dmfastqy`` (``-Y``)
single-end                 ``--fastqz`` (``-z``) ``--dmfastqz`` (``-Z``)
========================== ===================== =======================

For example, you could type ``-z {myfile}`` to input a FASTQ file of single-end
reads from a whole sample or ``-X {myfile1} -X {myfile2}`` to input two FASTQ
files of paired-end reads from a demultiplexed sample, replacing ``{myfile}``
with the actual path(s) of your file(s).

.. note::
    If you have paired-end reads in separate FASTQ files (``-x`` or ``-X``),
    then you must give both files of 1st and 2nd mates.
    See :ref:`fastq_pair` for instructions.

How to align one FASTQ file (single-end or interleaved paired-end reads)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

To align a FASTQ of single-end reads from a whole sample, use ``-z``::

    seismic align {refs.fa} -z {sample.fq.gz}

where ``{refs.fa}`` is the path to your FASTA file of reference sequences and
``{sample.fq.gz}`` is the path to your FASTQ file of the sample.

For a FASTQ of paired-end, interleaved reads that were demultiplexed, use ``-Y``
instead::

    seismic align {refs.fa} -Y {sample/ref.fq.gz}

where ``{sample/ref.fq.gz}`` is the path to your FASTQ file containing reads
from only one reference in the sample.

.. _fastq_pair:

How to align a pair of FASTQ files (paired-end reads in separate files)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

If your reads are paired-end and you have one FASTQ file containing all 1st
mates and another containing all 2nd mates, then you will need to provide both
FASTQ files.
There are two methods:

1.  Use ``-x``/``-X`` twice, once per FASTQ file::

        seismic align {refs.fa} -x {sample_R1.fq.gz} -x {sample_R2.fq.gz}

    where ``{sample_R1.fq.gz}`` and ``{sample_R2.fq.gz}`` are the paths to your
    FASTQ files of the 1st and 2nd mates, respectively.

2.  Make a new directory, move both FASTQ files into that directory, and provide
    the path to that directory with ``-x``/``-X``::

        mkdir {sample}
        mv {sample_R1.fq.gz} {sample_R2.fq.gz} {sample}
        seismic align {refs.fa} -x {sample}

    where ``{sample}`` is the new directory for both FASTQ files.

How to align multiple FASTQ files or pairs of paired-end FASTQ files
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

There are three ways to align multiple FASTQ files (or pairs thereof):

- **Use options more than once.**
  You can repeat any of ``-x``/``-y``/``-z``/``-X``/``-Y``/``-Z``, as well as
  mix them in one command.
  For example, to align one pair of paired-end FASTQ files (sample 1), one
  interleaved paired-end FASTQ file (sample 2), and two single-end FASTQ files
  (samples 3 and 4), you could type ::

    seismic align {refs.fa} -x {sample1_R1.fq.gz} -x {sample1_R2.fq.gz} -y {sample2.fq.gz} -z {sample3.fq.gz} -z {sample4.fq.gz}

  This method is most useful when you have a small number of FASTQ files.

- **Group FASTQ files of the same type into a directory.**
  Suppose you have 63 pairs of FASTQ files, with the files of mate 1s named
  ``sample-1_R1.fq.gz`` to ``sample-63_R1.fq.gz`` and the files of mate 2s named
  ``sample-1_R2.fq.gz`` to ``sample-63_R2.fq.gz``; plus demultiplexed single-end
  reads from three samples (I-III) and six references (A-F), named
  ``sample-I/ref-A.fq.gz`` to ``sample-III/ref-F.fq.gz``).
  You can align all of them with one command if you move the whole-sample,
  paired-end FASTQ files into their own directory, and the demultiplexed,
  single-end FASTQ files into another directory, and then give each directory
  after the appropriate options (``-x`` and ``-Z``, respectively)::

    mkdir {paired}
    mv sample-*_R?.fq.gz {paired}
    mkdir {dm-single}
    mv sample-I* {dm-single}
    seismic align {refs.fa} -x {paired} -Z {dm-single}

  This method is most useful when you have many FASTQ files.

- **Combine the first two methods.**
  Suppose you are working on two projects, have generated a set of many FASTQ
  files for each project, and want to process both sets.
  Currently, the FASTQ files for projects 1 and 2 are in directories ``proj1``
  and ``proj2``, and you want to keep them separate.
  You can process both directories with one command::

    seismic align {refs.fa} -x proj1 -x proj2

  This method is most useful when you have multiple directories of FASTQ files
  that you would like to keep separate.

.. note::
    If you give a directory for any of the FASTQ options, then SEISMIC-RNA will
    search for FASTQ files recursively, with no limit to the depth.

Align: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Align setting: Quality score encoding
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Your FASTQ files may encode quality scores in several schemes (for details, see
:ref:`phred_encodings`).
Modern Illumina sequencers use Phred+33 encoding, the default in SEISMIC-RNA.
To change the quality score encoding, use ``--phred-enc``.

.. note::
    If your FASTQ files do not use the Phred+33 encoding, then you must
    specify the correct Phred score encoding, or else Cutadapt and/or
    Bowtie 2 can produce incorrect output or fail outright.

If you do not know the encoding scheme of your FASTQ files, then you can process
them with `FastQC`_ and check the "Encoding" field in the "Basic Statisics" part
of the FastQC report:

- If the Encoding field says ``Illumina 1.0`` to ``1.7``, then your FASTQ files
  use Phred+64 encoding (``--phred-enc 64``).
- If the Encoding field says ``Illumina 1.8`` or greater, then your FASTQ files
  use Phred+33 encoding (``--phred-enc 33``, the default).
- Otherwise, you will need to search elsewhere for your encoding scheme.

Align setting: Quality assessment with FastQC
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To check the quality of your FASTQ files, SEISMIC-RNA runs `FastQC`_ by default.
To disable FastQC, use ``--no-fastqc``.
You can also enable automatic unzipping of the zipped output files from FastQC
with ``--qc-extract``.

Align setting: Trimming reads with Cutadapt
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To trim adapters and low-quality base calls before alignment, SEISMIC-RNA runs
`Cutadapt`_ by default.
To disable trimming, use ``--no-cut``.

How to trim adapter sequences
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Your reads may contain unwanted adapters (especially near their 3' ends), which
can cause problems such as misalignment (alignment to the wrong location).
Your adapter sequences depend on how your samples were prepared for sequencing
(i.e. on your library prep kit) and on your sequencing platform.
Since Illumina sequencers are the most widely used for mutational profiling,
SEISMIC-RNA defaults to the standard, minimal adapter sequences for Illumina
for both read 1 and (if paired-end) read 2:

- 5': ``GCTCTTCCGATCT``
- 3': ``AGATCGGAAGAGC``

If your samples have other adapters, then you can specify their sequences using

====== ====== ==============
 Side   Read   Option
====== ====== ==============
 5'     1      ``--cut-g1``
 5'     2      ``--cut-g2``
 3'     1      ``--cut-a1``
 3'     2      ``--cut-a2``
====== ====== ==============

.. _quality_trimming:

How to trim low-quality base calls
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Your reads may contain low-quality base calls (especially near their 3' ends),
which can cause misalignment and excessive mutations.
By default, SEISMIC-RNA trims base calls with quality scores less than 25, which
corresponds to a probability of 10\ :sup:`-2.5` = 0.3% that the base call is
incorrect (for an explanation, see :ref:`phred_encodings`).
You can set the quality threshold with ``--min-phred``.
We discourage using a quality threshold less than 25 because doing so could lead
to a background error rate that is too high for accurate mutational profiling
(e.g. 1% with ``--min-phred 20``), especially if you want to cluster your reads.

How to trim extra dark cycles (for Illumina two-channel chemistry)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Illumina sequencers using `two-channel chemistry`_ (e.g. NextSeq, NovaSeq, iSeq)
interpret the lack of color from either channel as G.
Consequently, if a DNA molecule is shorter than the read length, then the final
cycles of sequencing will produce no light and be `called as a string Gs`_.
Using ``--cut-nextseq`` tells Cutadapt to `trim high-quality Gs`_
from the 3' end of every read.
Trimming dark cycles can improve alignment in end-to-end mode, but it also trims
real G bases (which cannot be distinguished from artifactual ones) from the 3'
ends of reads.

How to further customize read trimming
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Refer to :ref:`cli_align` for the full list of options that SEISMIC-RNA can use
with Cutadapt, and the `Cutadapt reference guide`_ for details on each.
These options suffice for most users.
If you need more customization, then you can trim your FASTQ files externally
and then perform alignment within SEISMIC-RNA, using ``--no-cut`` to disable
additional trimming.

Align setting: Mapping reads with Bowtie 2
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA uses `Bowtie 2`_ to align your reads to your reference sequences.

How to pre-build a Bowtie 2 index (optional)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Bowtie 2 requires the FASTA file of reference sequences to be indexed.
You can have SEISMIC-RNA build the index for you automatically (the default) or
index your FASTA file manually.
With automatic indexing, SEISMIC-RNA builds the index in a temporary directory
and deletes it after alignment finishes.
This option is ideal for small sets of references (i.e. up to several hundred
sequences of several thousand nucleotides each) because building the index will
take on the order of seconds to minutes.
However, for large sets of references (e.g. an entire mammalian transcriptome),
building the index can take on the order of hours.
In this case, we recommend building the index yourself using the command ::

    bowtie2-build {refs}.fa {refs}

where ``{refs}.fa`` is the path of your FASTA file and ``{refs}`` is the path
without the FASTA file extension.
See the `Bowtie 2 Indexer manual`_ for more information on building an index.
Note that, while Bowtie 2 does not require the index to have the same name as
the FASTA file, SEISMIC-RNA does, so make sure that you use the same path for
the FASTA file and the index, except that the index path should not have the
FASTA file extension.

Indexing will generate six files with the extensions ``.1.bt2``, ``.2.bt2``,
``.3.bt2``, ``.4.bt2``, ``.rev.1.bt2``, and ``.rev.2.bt2``.
As long as all six files are in the same directory as and have the same name
(minus the file extension) as the FASTA file, SEISMIC-RNA will use the index.
Otherwise, SEISMIC-RNA will build and use a temporary index.

.. note::
    If you use a pre-built Bowtie 2 index, then SEISMIC-RNA does *not* verify
    that the index was actually built from the FASTA file of the same name.
    Discrepancies between the FASTA file and the index files can crash the Align
    and Relate steps or produce erroneous results.

How to choose between local and end-to-end alignment
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can align either whole reads (end-to-end mode) or the part of each read that
aligns best to the reference (local mode).
See the `description of alignment modes in Bowtie 2`_ for more details.

Generally, end-to-end mode yields spurious mutations (false positives) at the
ends of reads if the reads contain artifacts such as low-quality base calls or
untrimmed or improperly trimmed adapters.
Conversely, local mode misses real mutations (false negatives) within several
nucleotides of the ends of reads because such mutations cannot be part of the
best local alignment, which penalizes mutations and rewards matches.

For RNA mutational profiling, false positives generally cause more problems than
do false negatives, so SEISMIC-RNA uses local mode (``--bt2-local``) by default.
Use end-to-end mode (``--bt2-end-to-end``) only if you have a compelling reason
to do so (e.g. if you must quantify mutations at the ends of reads) and only
after carefully trimming any extraneous sequences from the ends of the reads.

How to align paired-end reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

If your reads are paired-end, then you have additional options for keeping or
discarding read pairs depending on how the two reads in the pair (called mates)
align relative to each other.
Bowtie 2 considers mates to align "concordantly" when their relative positions
match expectations and "discordantly" otherwise.
See the `Bowtie 2 manual for details on concordant/discordant alignments`_.
By default, SEISMIC-RNA keeps only concordantly aligned pairs.
To include discordantly aligned pairs too, add ``--bt2-discordant``.

Several options control which types of alignments are considered concordant
versus discordant.

You can specify where mates should align relative to each other: mates may
`overlap partially or completely, or dovetail`_.
By default, overlaps (partial and complete) are considered concordant, while
dovetailing is considered discordant.
You can treat complete overlaps as discordant with ``--bt2-no-contain``, or
dovetailed mates as concordant with ``--bt2-dovetail``.
Pairs that overlap partially (without dovetailing) are always concordant in
SEISMIC-RNA.

You can also specify the orientation of paired mates using ``--bt2-orient``.
The choices are ``fr`` (the 5'-most mate is forward, the 3'-most is reversed),
``rf`` (the 5'-most mate is reversed, the 3'-most is forward), or ``ff`` (both
mates are forward).
The default is ``fr`` (and if you are not sure which orientation you need, then
you probably need the default).

.. note::
    First, ``--bt2[-no]-contain``, ``--bt2[-no]-dovetail``, and ``--bt2-orient``
    choose which paired-end alignments count as concordant or discordant.
    If discordant, then ``--bt2-[no-]discordant`` choose whether to keep them.
    Using ``--bt2-no-contain`` and ``--bt2-discordant``, for example, would make
    alignments where one mate fully contains the other discordant (because of
    ``--bt2-no-contain``) but still kept (because of ``--bt2-discordant``),
    despite what the name "no-contain" would imply.

You can also enable `mixed mode`_ with ``--bt2-mixed``.
In mixed mode, if two mates fail to align as a pair, then Bowtie 2 will attempt
to align each mate individually, like a single-end read.
(It is possible in mixed mode for only one mate in a pair to align.)

How to filter aligned reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can filter alignments by `alignment score`_ and `mapping quality`_.

`Alignment score`_ measures how *well* a read aligns to *one specific location*
in *one reference sequence* and depends on the number of matches, substitutions,
and gaps, using the score parameters.
You can specify the minimum alignment score for local and end-to-end modes using
``--bt2-score-min-loc`` and ``--bt2-score-min-e2e``, respectively.
See the `section of the Bowtie 2 manual on alignment scores`_ for advice.

`Mapping quality`_ measures how *unique* an alignment is among *all locations*
in *all reference sequences*: high if the read aligns with a high alignment
score to exactly one location, low quality if it aligns with similar alignment
scores to multiple locations in the reference (and thus it is hard to determine
a single location where the read aligns).
The default minimum mapping quality is 25, meaning that the probability that the
chosen location is incorrect is 10\ :sup:`-2.5` = 0.3%.
You can change the minimum mapping quality using ``--min-mapq``.

How to filter by number of aligned reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Alignment maps containing very few reads are not generally useful for mutational
profiling, due to low coverage per position.
When aligning to many references (e.g. an entire transcriptome), most references
will receive few reads, producing many output files that would be unusable for
further processing.
To prevent unusable files from cluttering your output directory, you can choose
to have alignment map files with insufficient reads deleted automatically.
The default minimum is 1000 reads, which you can change using ``--min-reads``.
With no minimum (``--min-reads 0``), no files are deleted automatically.

How to further customize alignment
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

See :ref:`cli_align` for the full list of options that SEISMIC-RNA can use with
Bowtie 2, and the `Bowtie 2 manual`_ for details on each of these options.
These options suffice for most users.
If you need more customization, then you can align your FASTQ files externally
and pass the alignment maps into SEISMIC-RNA during :doc:`./relate`.

.. _bam_vs_cram:

Align setting: Format of alignment maps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can choose to output alignment map files in either BAM or CRAM format.
For information on these file formats, see :doc:`../../formats/data/xam`.
The default is CRAM (``--cram``); you can switch to BAM using ``--bam``.

Alignment maps in CRAM format are smaller than their BAM counterparts, and hence
better suited to long-term storage.
However, the better compression of CRAM files comes at three costs:

- A CRAM file must be accompanied by a FASTA file storing the sequence of every
  reference that appears in the header of the CRAM file.
  A CRAM file stores only the relative path to its FASTA file, not the sequence
  information, which enables the CRAM file to be much smaller than it would be
  if it did need to contain its own sequences.
  Because the FASTA file existed before and during the alignment, having this
  FASTA file accompany the CRAM file usually incurs no extra cost.
  However, moving or deleting the FASTA will break the CRAM file.
  As a safeguard against this fragility, SEISMIC-RNA keeps a copy of the FASTA
  file in the same directory as the output CRAM file.
  Creating an actual copy would require more storage space and defeat the point
  of CRAM's smaller file size, so SEISMIC-RNA actually makes a `hard link`_ --
  not a copy -- which requires minimal extra space.
  In some circumstances, making a hard link can fail, in which case SEISMIC-RNA
  will resort to copying the FASTA file instead.
- Reading and writing CRAM files is slower than for BAM files due to the extra
  effort needed for compressing and decompressing CRAM files.
- In the `CIGAR strings`_, distinction between reference matches (``=``) and
  substitutions (``X``) is lost upon compressing to CRAM format.
  Thus, the Relate step must perform extra work to determine if each non-indel
  position is a match or substitution, which makes it run more slowly than it
  would if the distinction had been preserved.

In general, use CRAM format if minimizing the size of your alignment map files
is a priority, especially for long-term storage.
Use BAM format to make the ``align`` and ``relate`` steps run faster, and to
make the output files more portable (since BAM files are self-contained, while
CRAM files will break without the FASTA file that accompanies them).

Align: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files except FastQC reports are written to ``{out}/{sample}/align``,
where ``{out}`` is your output directory and ``{sample}`` is the sample name.

Align output file: FastQC reports
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If you run FastQC, then its report files go into ``{out}/{sample}/qc``.
The directory ``{out}/{sample}/qc/initial`` contains the FastQC reports for your
initial FASTQ files, before trimming.
If you also run trimming, then reports for the post-trimmed FASTQ files go into
``{out}/{sample}/qc/trimmed``.

In each directory, FastQC outputs ``{fq_name}_fastqc.html`` (the FastQC report)
and ``{fq_name}_fastqc.zip`` (extra information), where ``{fq_name}`` comes from
the original FASTQ file.
If you add ``--qc-extract``, then each ``{fq_name}_fastqc.zip`` will be unzipped
to the directory ``{fq_name}_fastqc``.
For details on these outputs, see the documentation for `FastQC`_.

Align output file: Alignment maps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Alignment maps store the location in the reference sequence to which each read
aligned, plus the Phred quality scores, mapping quality, and mutated positions.
(For more information on alignment maps, see :doc:`../../formats/data/xam`.)
SEISMIC-RNA outputs alignment maps where every read aligns to the same reference
(although this is not a restriction outside of SEISMIC-RNA).
Each alignment map is written to ``{ref}.{xam}``, where ``{ref}`` is the name of
the reference to which the reads aligned, and ``{xam}`` is the file extension
(depending on the selected format).
SEISMIC-RNA can output alignment maps in either BAM or CRAM format.
For a comparison of these formats, see :ref:`bam_vs_cram`.

Align output file: Reference sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If you choose to output alignment maps in CRAM format, then you also get a FASTA
file(s) of the reference sequence(s) alongside the CRAM files.
If the reads came from a whole sample, then a single FASTA file with the same
name as the input FASTA file will be output.
The output file will be a `hard link`_ to the input file, if possible, to avoid
consuming unnecessary storage space.
If the reads were demultiplexed before alignment, then for each CRAM file, a
FASTA file with the same name (up to the file extension) will be output.
In both cases, each FASTA file will be indexed using `samtools faidx`_ to speed
up reading the CRAM files.
If you choose to output alignment maps in BAM format, then you get (and need)
no FASTA files alongside them.

.. _align_unaligned:

Align output file: Unaligned reads
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In addition to the alignment maps, SEISMIC-RNA outputs FASTQ file(s) of reads
that Bowtie 2 could not align:

- Each whole-sample FASTQ file of single-end (``-z``) or interleaved (``-y``)
  reads yields one file: ``unaligned.fq.gz``
- Each pair of whole-sample FASTQ files of 1st and 2nd mates (``-x``) yields two
  files: ``unaligned.fq.1.gz`` and ``unaligned.fq.2.gz``
- Each demultiplexed FASTQ file of single-end (``-Z``) or interleaved (``-Y``)
  reads yields one file: ``{ref}__unaligned.fq.gz``
- Each pair of demultiplexed FASTQ files of 1st and 2nd mates (``-X``) yields
  two files: ``{ref}__unaligned.fq.1.gz`` and ``{ref}__unaligned.fq.2.gz``

where ``{ref}`` is the reference for demultiplexed FASTQ files.

You can disable outputting unaligned using ``--bt2-no-un``.

Align output file: Align report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA also writes a report file, ``align-report.json``, that records the
settings you used for running the Align step and summarizes the results.
See :doc:`../../formats/report/align` for more information.

Check the number of reads that aligned overall
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Under "Number of reads after alignment", the report says how many single-end
and/or paired-end reads were in the FASTQ file(s), and how many reads aligned.
This information is copied verbatim from the `alignment summary`_ of Bowtie 2;
see its documentation for more details.
For paired-end reads, each pair counts as one read.

Check the number of reads that aligned to each reference
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Under "Number of reads aligned by reference", the report lists every reference
in your input FASTA file and the number of reads that aligned to it.
For paired-end reads, each pair counts as one read.

Align: Troubleshoot and optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Align produces alignment map files too slowly
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

First, try running the Align step using more processors (with ``--max-procs``),
at the cost of using more memory.
If, as a result, :ref:`align_crash_hang`, then try adjusting the settings of
Bowtie 2 to increase the speed, at the risk of overlooking valid alignments.
See :ref:`cli_align` for the Bowtie 2 settings you can adjust in SEISMIC-RNA,
and the `Bowtie 2 manual`_ for more detailed descriptions.

.. _align_crash_hang:

Align crashes or hangs without producing alignment map files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Most likely, your system has run out of memory.
You can confirm using a program that monitors memory usage (such as ``top`` in a
Linux/macOS terminal, Activity Monitor on macOS, or Task Manager on Windows).
If so, then rerun the Align step using fewer processors (with ``--max-procs``)
to limit the memory usage, at the cost of slower alignment.

Fewer reads aligned than you expected
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Try the following steps (in this order):

1.  Ensure you are using Bowtie version 2.5.1 or later (version 2.5.0 has a bug
    that affects alignment rate).
    You can check the version with ``bowtie2 --version | head -n 1``.
2.  Double check that your FASTA file has the correct reference sequence(s) and
    that, if you pre-built the Bowtie 2 index before running ``seismic align``,
    that you indexed the correct FASTA file.
3.  Examine the reads that failed to align (see :ref:`align_unaligned`).
    Choose several reads randomly, copy one or two 20 - 40 nt segments from the
    middle of each read, and check if the segments come from any known sources
    by querying `BLAST`_ (or similar tools).
    Identifying the sources of unaligned reads can help determine the cause of
    the problem (e.g. contamination with ribosomal RNA or foreign nucleic acids
    such as from *Mycoplasma*) and whether the reads that did align are usable.

.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Cutadapt: https://cutadapt.readthedocs.io/en/stable/
.. _two-channel chemistry: https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology/2-channel-sbs.html
.. _called as a string Gs: https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/
.. _trim high-quality Gs: https://cutadapt.readthedocs.io/en/stable/guide.html#nextseq-trim
.. _Cutadapt reference guide: https://cutadapt.readthedocs.io/en/stable/reference.html
.. _Bowtie 2: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _Bowtie 2 Indexer manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
.. _description of alignment modes in Bowtie 2: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment
.. _alignment score: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar
.. _section of the Bowtie 2 manual on alignment scores: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#valid-alignments-meet-or-exceed-the-minimum-score-threshold
.. _mapping quality: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mapping-quality-higher-more-unique
.. _CIGAR strings: https://samtools.github.io/hts-specs/
.. _view command in Samtools: https://www.htslib.org/doc/samtools-view.html
.. _Bowtie 2 manual for details on concordant/discordant alignments: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#concordant-pairs-match-pair-expectations-discordant-pairs-dont
.. _alignment summary: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#alignment-summary
.. _mixed mode: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mixed-mode-paired-where-possible-unpaired-otherwise
.. _overlap partially or completely, or dovetail: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other
.. _Bowtie 2 manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
.. _hard link: https://en.wikipedia.org/wiki/Hard_link
.. _samtools faidx: https://www.htslib.org/doc/samtools-faidx.html
