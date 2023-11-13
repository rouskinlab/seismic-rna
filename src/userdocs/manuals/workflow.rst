
Running the Workflow
========================================================================

There are two points from which you can begin the main workflow:

- If you are starting from files of raw sequencing reads (FASTQ format),
  then begin at the step :ref:`wf_align`.
- If you are starting from files of aligned sequencing reads (SAM, BAM,
  or CRAM format), then begin at the step :ref:`wf_relate`.

.. note::
    The command ``seismic all`` accepts both types of inputs and runs
    the entire workflow.
    See :ref:`wf_all` for more details.


.. _wf_align:

Align sequencing reads to reference sequences
------------------------------------------------------------------------

Align: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Align input file: Reference sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Alignment requires exactly one file of reference sequences, which must
be DNA sequences (only A, C, G, T, and N) in FASTA format.
If needed, you may clean the FASTA file with the :doc:`./faclean` tool.
See :doc:`../formats/data/fasta` for details.

Align input file: Sequencing reads
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Alignment requires one or more files of sequencing reads, which must be
DNA sequences (only A, C, G, T, and N) in FASTQ format.
See :doc:`../formats/data/fastq` for details.

Sequencing reads can be single-end or paired-end
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Sequencing reads can be single-end or paired-end (for more details, see
:ref:`fastq_endedness`).
The SEISMIC-RNA workflow can handle both single- and paired-end reads,
provided that no FASTQ file contains a mixture of both types of reads.

Sequencing reads can come from whole or demultiplexed samples
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

FASTQ files can contain reads from a whole sample (which can possibly
come from many different reference sequences) or can be demultiplexed
before alignment so that they contain reads from one reference sequence.

How to specify the endedness and source of sequencing reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Specify the endedness (single-end, paired-end interleaved, or paired-end
separated) and source of the reads (whole or demultiplexed sample) in a
FASTQ file by giving the file name using the appropriate option below.
Each option has a "long form" (starting with ``--``) and a "short form"
(starting with ``-`` and comprising a single character), which are both
indicated in this table, separated by a ``/``:

========================== =============== =================
Endedness                  Whole Sample    Demultiplexed
========================== =============== =================
paired-end, separate files ``--fastqx/-x`` ``--dmfastqx/-X``
paired-end, interleaved    ``--fastqy/-y`` ``--dmfastqy/-Y``
single-end                 ``--fastqz/-z`` ``--dmfastqz/-Z``
========================== =============== =================

.. note::
    As with all options, the long and short forms are interchangable.
    For brevity, the rest of this manual uses only the short forms.

.. note::
    If using paired-end reads in separate FASTQ files, both files of 1st
    and 2nd mates are required. See :ref:`fastq_pair` for instructions.

How to align one FASTQ file (single-end or interleaved paired-end reads)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

To align a FASTQ of single-end reads from a whole sample, use ``-z``::

    seismic align {refs.fa} -z {sample.fq.gz}

where ``{refs.fa}`` is the path to the FASTA file of reference sequences
and ``{sample.fq.gz}`` is the path to the FASTQ file of the sample.

For a FASTQ of paired-end, interleaved reads that were demultiplexed,
use ``-Y`` instead::

    seismic align {refs.fa} -Y {sample/ref.fq.gz}

where ``{sample/ref.fq.gz}`` is the path to the FASTQ file containing
reads from only one reference in the sample.

.. _fastq_pair:

How to align a pair of FASTQ files (paired-end reads in separate files)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

If your reads are paired-end and you have one FASTQ file containing all
1st mates and another containing all 2nd mates, then you will need to
provide both FASTQ files.
There are two methods:

1.  Use the option ``-x`` (or ``-X``) twice, once per FASTQ file::

        seismic align {refs.fa} -x {sample_R1.fq.gz} -x {sample_R2.fq.gz}

    where ``{sample_R1.fq.gz}`` and ``{sample_R2.fq.gz}`` are the paths
    to the FASTQ files of the 1st and 2nd mates, respectively.

2.  Make a new directory, move both FASTQ files into that directory, and
    provide the path to that directory with ``-x`` (or ``-X``)::

        mkdir {sample}
        mv {sample_R1.fq.gz} {sample_R2.fq.gz} {sample}
        seismic align {refs.fa} -x {sample}

    where ``{sample}`` is the new directory for both FASTQ files.

How to align multiple FASTQ files or pairs of paired-end FASTQ files
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

There are three ways to align multiple FASTQ files (or pairs thereof):

1.  Use options more than once.
    The options for FASTQ files can all be given multiple times, and can
    even be mixed in one command.
    For example, to align a pair of paired-end FASTQ files (sample 1),
    an interleaved paired-end FASTQ file (sample 2), and two single-end
    FASTQ files (samples 3 and 4), use the following options::

        seismic align {refs.fa} -x {sample1_R1.fq.gz} -x {sample1_R2.fq.gz} -y {sample2.fq.gz} -z {sample3.fq.gz} -z {sample4.fq.gz}

    This method is most useful when you have a few FASTQ files.

2.  Group FASTQ files of the same type into a directory.
    For example, suppose you have 63 FASTQ files each of paired-end 1st
    mates (named ``sample-1_R1.fq.gz`` to ``sample-63_R1.fq.gz``) and
    2nd mates (named analogously but with ``R2``), plus demultiplexed
    single-end reads from three samples (I-III) and six references (A-F)
    (named ``sample-I/ref-A.fq.gz`` to ``sample-III/ref-F.fq.gz``).
    Move the separate paired-end FASTQ files into their own directory,
    and the demultiplexed single-end FASTQ files into another directory,
    then run alignment by passing each directory of FASTQ files::

        mkdir {paired}
        mv sample-*_R?.fq.gz {paired}
        mkdir {dm-single}
        mv sample-I* {dm-single}
        seismic align {refs.fa} -x {paired} -Z {dm-single}

    This method is most useful when you have many FASTQ files.

3.  Combine methods 1 and 2.
    Suppose you are working on two projects, have generated a set of
    many FASTQ files for each project, and want to process both sets.
    Currently, the FASTQ files for projects 1 and 2 are in directories
    ``proj1`` and ``proj2``, and you want to keep them separate.
    You can process both directories with one command::

        seismic align {refs.fa} -x proj1 -x proj2

    This method is most useful when you have multiple directories of
    FASTQ files that you would like to keep separate.

.. note::
    If a directory is given for any of the FASTQ options, then it will
    be searched for FASTQ files recursively, with no limit to the depth.
    Thus, the given directory can have deeply nested subdirectories, and
    SEISMIC-RNA will still find and process any FASTQ files within them.

Align: Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Align option: Phred score encoding
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA defaults to using Phred+33 encoding for FASTQ files, which
is standard on modern Illumina sequencers.
To change the Phred score encoding, use the option ``--phred-enc``.
See :ref:`phred_encodings` for more information.

.. note::
    If your FASTQ files do not use the Phred+33 encoding, then you must
    specify the correct Phred score encoding, or else Cutadapt and/or
    Bowtie 2 can produce incorrect output or fail outright.

If you do not know the encoding scheme of your FASTQ files, then you may
be able to determine it by using `FastQC`_ or ``seismic align`` (which
runs FastQC automatically).
In the HTML report file generated by FastQC, check the "Encoding" field
in the "Basic Statisics" section:

- If the Encoding field says ``Illumina 1.0`` to ``1.7``, then your
  FASTQ files use Phred+64 encoding (``--phred-enc 64``).
- If the Encoding field says ``Illumina 1.8`` or greater, then your
  FASTQ files use Phred+33 encoding (``--phred-enc 33``, the default).
- Otherwise, you will need to search elsewhere for your encoding scheme
  to determine the Phred score offset.

Align option: Quality assessment with FastQC
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, each FASTQ file is processed with `FastQC`_, both before and
after trimming, in order to find any potential problems.
FastQC can be disabled with the flag ``--no-fastqc``.
To enable automatic extraction of the zipped output files from FastQC,
add the flag ``--qc-extract``.

Align option: Trimming reads with Cutadapt
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, each FASTQ file and pair of mated FASTQ files is trimmed for
adapters and low-quality bases using `Cutadapt`_. To disable trimming,
add the flag ``--no-cut``.

How to trim adapter sequences
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

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

How to trim low-quality base calls
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Base calls on either end of a read that fall below a minimum Phred score
quality are trimmed with Cutadapt.
The default minimum quality is 25, which corresponds to a probability of
1 - 10 :sup:`-2.5` = 0.997 that the base call is correct.
(See :ref:`phred_encodings` for more details).
To change the quality threshold, use the option ``--min-phred``.

How to use Cutadapt to trim dark cycles (for Illumina NextSeq or iSeq)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

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

How to further customize read trimming
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

See :ref:`cli_align` for the full list of options that SEISMIC-RNA can
use with Cutadapt, and the `Cutadapt reference guide`_ for details on
each of these options.
These options should suffice for most users.
If you require a more customized adapter trimming workflow, you can trim
your FASTQ files outside of SEISMIC-RNA, then perform alignment within
SEISMIC-RNA, using the option ``--no-cut`` to disable additional adapter
trimming.

Align option: Mapping reads with Bowtie 2
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

How to pre-build a Bowtie 2 index (optional)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Bowtie 2 requires the FASTA file of reference sequences to be indexed.
A Bowtie 2 index comprises six files, all in the same directory as and
with the same name as the FASTA file, with the extensions ``.1.bt2``,
``.2.bt2``, ``.3.bt2``, ``.4.bt2``, ``.rev.1.bt2``, and ``.rev.2.bt2``.

If the index is missing, then SEISMIC-RNA will create a temporary index
automatically using ``bowtie2-build`` each time you run alignment.
Automatic indexing is efficient when the FASTA file is small: several
hundred reference sequences or fewer.
For larger FASTA files (e.g. a whole eukaryotic transcriptome), building
a temporary index each time alignment is run becomes costly.
In this case, it is more efficient to pre-build the index, which you can
do with this command::

    bowtie2-build {refs}.fa {refs}

where ``{refs}`` is the path to and name of your FASTA file.
See the `Bowtie 2 Indexer manual`_ for more details.

.. note::
    If you use a pre-built Bowtie 2 index, then SEISMIC-RNA does *not*
    verify that the index was actually built from the FASTA file of the
    same name.
    You can assume the index is correct if you build it using the above
    command and avoid modifying or replacing the FASTA and index files.
    Discrepancies between the FASTA file and the index files can crash
    the ``align`` and ``relate`` steps or produce erroneous results.

How to choose between local and end-to-end alignment
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

During alignment, Bowtie 2 can either align the entirety of each read
(end-to-end mode) or find and align only the section of the read that
yields the best alignment score (local mode).
See the `description of alignment modes in Bowtie 2`_ for more details.

Generally, end-to-end mode yields spurious mutations (false positives)
at the ends of reads if the reads contain artifacts such as low-quality
base calls or untrimmed or improperly trimmed adapters.
Conversely, local mode misses real mutations (false negatives) within
several nucleotides of the ends of reads because such mutations are not,
by definition, part of the best local alignment.

Concerning RNA mutational profiling, false positives are generally much
more problematic than false negatives, so SEISMIC-RNA uses local mode
(``--bt2-local``) by default.
Use end-to-end mode (``--bt2-end-to-end``) only if you have a compelling
reason to do so (e.g. if it is essential to detect mutations at the ends
of reads) and only after carefully trimming any extraneous sequences
from the ends of the reads.

How to align paired-end reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

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

How to filter aligned reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Alignments can be filtered by `alignment score`_ and `mapping quality`_,
which are distinct properties.

`Alignment score`_ measures how well a read aligns to a given location
in the reference.
It is calculated from the number of matches, substitutions, and gaps
using the score parameters.
The minimum alignment scores for local and end-to-end modes can be set
using ``--bt2-score-min-loc`` and ``--bt2-score-min-e2e``, respectively.
See the `section of the Bowtie 2 manual on alignment scores`_ for advice
on setting this parameter.

`Mapping quality`_ measures how unique an alignment is: high quality if
the read aligns with a high score to exactly one location, low quality
if it aligns with similar scores to multiple locations in the reference.
The default minimum quality is 25, which corresponds to a confidence of
1 - 10 :sup:`-2.5` = 0.997 that the read has aligned correctly.
To change the quality threshold, use the option ``--min-mapq``.
For those searching for this option in Bowtie 2, you will not find it.
Instead, reads with insufficient mapping quality are filtered out after
alignment using the `view command in Samtools`_.

How to filter by number of aligned reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

In general, alignment maps containing very few reads are not useful for
mutational profiling, due to their inherently low coverage per position.
Worse, if aligning to a very large number of references (e.g. an entire
transcriptome), most of the references would likely receive insufficient
reads, so most of the (many) output XAM files would be useless clutter.

To remedy this inconvenience, after alignment has finished, XAM files
with fewer than a minimum number of reads are automatically deleted.
The default is 1000, which can be set using the option ``--min-reads``.
Setting ``--min-reads`` to 0 disables automatically deleting XAM files.

How to further customize alignment
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

See :ref:`cli_align` for the full list of options that SEISMIC-RNA can
use with Bowtie 2, and the `Bowtie 2 manual`_ for details on each of
these options.
These options should suffice for most users.
If you require a more customized alignment workflow, you can align your
your FASTQ files outside of SEISMIC-RNA, then pass the resulting XAM
files into SEISMIC-RNA at the step :ref:`wf_relate`.


.. _bam_vs_cram:

Align option: Format of alignment maps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA can output alignment map files in either BAM or CRAM format.
For details on these file formats, see :doc:`../../formats/data/xam`.
The default is CRAM format (option ``--cram``); BAM format is enabled
using the option ``--bam``.

Alignment maps in CRAM format are smaller than their BAM counterparts,
and hence better suited to long-term storage.
However, the better compression of CRAM files comes at three costs:

- A CRAM file must be accompanied by a FASTA file storing the sequence
  of every reference that appears in the header of the CRAM file.
  A CRAM file stores only the relative path to its FASTA file, not the
  sequence information, which enables the CRAM file to be much smaller
  than it would be if it did need to contain its own sequences.
  Because the FASTA file existed before and during the alignment, having
  this FASTA file accompany the CRAM file usually incurs no extra cost.
  However, moving or deleting the FASTA will break the CRAM file.
  As a safeguard against this fragility, SEISMIC-RNA keeps a copy of the
  original FASTA file in the same directory as the output CRAM file.
  Creating an actual copy would require more storage space and defeat
  the purpose of CRAM's smaller file size, so SEISMIC-RNA actually makes
  a `hard link`_ -- not a copy -- which requires minimal extra space.
  In some circumstances, making a hard link can fail, in which case
  SEISMIC-RNA will resort to copying the FASTA file instead.
- Reading and writing CRAM files is slower than for BAM files due to the
  extra effort needed for compressing and decompressing CRAM files.
- In the `CIGAR strings`_, distinction between reference matches (``=``)
  and substitutions (``X``) is lost upon compressing to CRAM format.
  Thus, ``seismic relate`` must perform extra work to determine if each
  non-gapped position is a match or substitution, which makes it run
  more slowly than it would if the distinction had been preserved.

In general, use CRAM format if minimizing the size of your alignment
map files is a priority, especially for long-term storage.
Use BAM format to make the ``align`` and ``relate`` steps run faster,
and to increase the robustness of the output files (because BAM files
are self-contained, while CRAM files will break without the FASTA file
that accompanies them).

Align: output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Align output file: FastQC reports
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If FastQC is run, then it outputs files to ``{out}/{sample}/qc``, where
``{out}`` is the output directory (``--out-dir``) and ``{sample}`` is
the name of the sample.
The directory ``{out}/{sample}/qc/init`` is always created and contains
FastQC reports of the initial FASTQ files.
If adapter/quality trimming was run, ``{out}/{sample}/qc/trim`` is also
created for FastQC reports of the trimmed FASTQ files.

In each directory (``init`` and ``trim``), FastQC writes two files for
each FASTQ file: ``{fq_name}_fastqc.html`` and ``{fq_name}_fastqc.zip``,
where ``{fq_name}`` is the name of the original FASTQ file up to the
file extension.
If the option ``--qc-extract`` is given, then FastQC will also unzip
``{fq_name}_fastqc.zip`` to the directory ``{fq_name}_fastqc``.
For details on these outputs, see the documentation for `FastQC`_.

Align output file: Alignment maps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The most important outputs of ``seismic align`` are alignment map files.
Alignment maps store the location in the reference sequence to which
each read aligned, as well as the Phred quality scores, mapping quality,
and mutated positions.
SEISMIC-RNA outputs alignment maps where every read aligns to the same
reference (although this is not a restriction outside of SEISMIC-RNA).
Each alignment map is written to ``{out}/{sample}/align/{ref}.{xam}``,
where ``{out}`` is the output directory (``--out-dir``), ``{sample}`` is
the name of the sample from which the reads came, ``{ref}`` is the name
of the reference to which the reads aligned, and ``{xam}`` is the file
extension (depending on the selected format).
SEISMIC-RNA can output alignment maps in either BAM or CRAM format.
For a comparison of these formats, see :ref:`bam_vs_cram`.

Align output file: Reference sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If the alignment maps are output in CRAM format, then FASTA file(s) of
the reference sequence(s) are also output alongside the CRAM files.
If the sequencing reads came from a whole sample, then a single FASTA
file, bearing the same name as the input FASTA file, will be output.
The output file will be a `hard link`_ to the input file, if possible,
to avoid consuming unnecessary storage space.
If the sequencing reads were demultiplexed before alignment, then for
each output CRAM file, a FASTA file with the same name (up to the file
extension) will be written to the same directory.
In both cases, each output FASTA will be indexed using `samtools faidx`_
to speed up reading the CRAM files.
If the alignment maps are output in BAM format, then FASTA files are not
output alongside them.

.. _wf_unaligned:

Align output file: Unaligned reads
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In addition to the alignment maps, SEISMIC-RNA outputs FASTQ file(s) of
reads that Bowtie 2 could not align to ``{out}/{sample}/align``:

- Each whole-sample FASTQ file of single-end (``-z``) or interleaved
  (``-y``) reads yields one file: ``unaligned.fq.gz``
- Each pair of whole-sample FASTQ files of 1st and 2nd mates (``-x``)
  yields two files: ``unaligned.fq.1.gz`` and ``unaligned.fq.2.gz``
- Each demultiplexed FASTQ file of single-end (``-Z``) or interleaved
  (``-Y``) reads yields one file: ``{ref}__unaligned.fq.gz``
- Each pair of demultiplexed FASTQ files of 1st and 2nd mates (``-X``)
  yields two files:
  ``{ref}__unaligned.fq.1.gz`` and ``{ref}__unaligned.fq.2.gz``

where ``{ref}`` is the reference for demultiplexed FASTQ files.

Outputting these files of unaligned reads can be disabled using the
option ``--bt2-no-un``.

Align output file: Report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A report file is written that records the settings used to run alignment
and summarizes the results of alignment.
See :doc:`../formats/report/align` for more information.

Align: Troubleshooting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Troubleshooting a lower-than-expected alignment rate
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If the percent of reads aligning to the reference is less than expected,
then try the following steps (in this order):

1.  Ensure you are using Bowtie version 2.5.1 or later (version 2.5.0
    has a bug that affects alignment rate).
    You can check the version with ``bowtie2 --version | head -n 1``.
2.  Double check that the FASTA has the correct reference sequence(s)
    and that, if the Bowtie 2 index was pre-built before the align step,
    that the correct FASTA file was used.
3.  Examine the reads that failed to align (see :ref:`wf_unaligned`).
    Choose several reads randomly and check if they could come from any
    known sources by querying `BLAST`_ (or similar tools) for short
    (20 - 40 nt) segments of each read.
    Identifying the sources of unaligned reads can help determine the
    cause of the problem (e.g. contamination with ribosomal or foreign
    RNA such as from *Mycoplasma*, incorrect indexes used during FASTQ
    generation) and whether the reads that did align are still usable.

.. _wf_relate:

Relate each read to every reference position
------------------------------------------------------------------------

Relate: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Relate input file: Reference sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Relate requires exactly one file of reference sequences, which must be
DNA sequences (only A, C, G, T, and N) in FASTA format.
If needed, you may clean the FASTA file with the :doc:`./faclean` tool.
See :doc:`../formats/data/fasta` for details.

Relate input file: Alignment maps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Relate requires one or more alignment map files, each of which must be
in SAM, BAM, or CRAM format (collectively, "XAM" format).
See :doc:`../formats/data/xam` for details.

.. note::
    The references in the FASTA file must match those to which the reads
    in the alignment map were aligned.
    Discrepancies can cause the ``relate`` command to fail or produce
    erroneous relation vectors.
    This problem will not occur if you use the same (unaltered) FASTA
    file for both the ``align`` and ``relate`` commands, or run both
    at once using the command ``seismic all``.

List every alignment map file after the FASTA file.
Refer to :doc:`./inputs` for details on how to list multiple files.
For example, to compute relation vectors for reads from ``sample-1``
aligned to references ``ref-1`` and ``ref-2``, and from ``sample-2``
aligned to reference ``ref-1``, use the following command::

    seismic relate {refs.fa} sample-1/align/ref-1.cram sample-1/align/ref-2.cram sample-2/align/ref-1.cram

where ``{refs.fa}`` is the path to the file of reference sequences.

Relate: Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Relate options shared with alignment
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Because this workflow can be started from the ``align`` or ``relate``
commands, the latter duplicates some of the options of the former:
``--phred-enc``, ``--min-mapq``, ``--min-reads``, and ``--out-dir`` have
the same functions in ``relate`` and ``align``.

Relate option: Minimum Phred score
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Like ``align``, ``relate`` also has the option ``--min-phred``, but its
meaning is different than that during the ``align`` step.
In ``relate``, base calls with Phred scores below ``--min-phred`` are
considered ambiguous matches or substitutions, as if they were ``N``s.
For example, if the minimum Phred score is 25 (the default) and a base
``T`` is called as a match with a Phred score of 20, then it would be
marked as possibly a match and possibly a subsitution to A, C, or G.
See :doc:`../data/relate` for more information.

Relate option: Ambiguous insertions and deletions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The most tricky problem in computing relation vectors is that insertions
and deletions ("indels") in repetitive regions cause ambiguities.
SEISMIC-RNA introduces a new algorithm for identifying ambiguous indels
(see :doc:`../algos/ambrel` for more information).
This algorithm is enabled by default.
If it is not necessary to identify ambiguous indels, then the algorithm
can be disabled with ``--no-ambrel``, which will speed up ``relate`` at
the cost of reducing its accuracy on indels.

Relate option: Batch size
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For an explanation of batching and how to use it, see :ref:`batches`.

The dataset is partitioned into batches by the ``relate`` command.
The option ``--batch-size`` sets a target amount of data for each batch,
in millions of base calls (megabases).
This calculation considers the total number of relationships per read,
which equals the length of the reference sequence.
Thus, the number of base calls *B* is the product of the number of reads
*N* and the length of the reference sequence *L*:

*B* = *NL*

Since *L* is known and ``--batch-size`` specifies a target size for *B*,
*N* can be solved for:

*N* = *B*/*L*

SEISMIC-RNA will aim to put exactly *N* reads in each batch but the last
(the last batch can be smaller because it has just the leftover reads).
If the reads are single-ended or if alignment was not run in mixed mode,
then every batch but the last will contain exactly *N* reads.
If mixed mode was used, then batches may contain more than *N* reads, up
to a maximum of 2 *N* in the extreme case that every read in the batch
belonged to a pair in which the other mate did not align.

Relate: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Relate output file: Batch of relation vectors
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



Relate output file: Report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A report file is written that records the settings used to generate
relation vectors summarizes the results.
See :doc:`../formats/report/relate` for more information.

.. _wf_all:

Run the entire workflow with one command
------------------------------------------------------------------------

.. note::
    ``seismic all`` accepts FASTQ, SAM/BAM/CRAM, relate/mask/cluster report, and
    table files and directories as inputs.

From BAM, report, and/or table file(s)::

    seismic all refs.fa out/sample/align/Ref.bam out/sample/*/*-report.json out/sample/table/*/*.csv


.. note::
    Only the align, relate, mask, and table steps run by default. Enable
    clustering by specifying ``--max-clusters`` (``-k``) followed by the
    maximum number of clusters to attempt. Enable structure prediction
    with the flag ``--fold``.

.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Cutadapt: https://cutadapt.readthedocs.io/en/stable/
.. _Cutadapt reference guide: https://cutadapt.readthedocs.io/en/stable/reference.html
.. _Bowtie 2 Indexer manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
.. _description of alignment modes in Bowtie 2: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment
.. _alignment score: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar
.. _section of the Bowtie 2 manual on alignment scores: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#valid-alignments-meet-or-exceed-the-minimum-score-threshold
.. _mapping quality: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mapping-quality-higher-more-unique
.. _CIGAR strings: https://samtools.github.io/hts-specs/
.. _view command in Samtools: https://www.htslib.org/doc/samtools-view.html
.. _Bowtie 2 manual for details on concordant/discordant alignments: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#concordant-pairs-match-pair-expectations-discordant-pairs-dont
.. _mixed mode: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mixed-mode-paired-where-possible-unpaired-otherwise
.. _overlap partially or completely, or dovetail: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other
.. _Bowtie 2 manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
.. _hard link: https://en.wikipedia.org/wiki/Hard_link
.. _samtools faidx: https://www.htslib.org/doc/samtools-faidx.html
.. _glob patterns: https://en.wikipedia.org/wiki/Glob_(programming)
