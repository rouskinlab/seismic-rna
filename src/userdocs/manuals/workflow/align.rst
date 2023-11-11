
Align the sequencing reads
------------------------------------------------------------------------

Input file: FASTA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alignment requires exactly one file of reference sequences, which must
be in FASTA format.
See :doc:`../../formats/data/fasta` for details.
You can use the  need to clean up a FASTA file that


Input file(s): FASTQ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alignment requires one or more files of sequencing reads, which must be
in FASTQ format.
See :doc:`../../formats/data/fastq` for details.


FASTQ files of single-end versus paired-end reads
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Sequencing reads can be single-end or paired-end (for more details, see
:ref:`fastq_endedness`).
The SEISMIC-RNA workflow can handle both single- and paired-end reads,
provided that no FASTQ file contains a mixture of both types of reads.


FASTQ files of whole samples versus demultiplexed samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

FASTQ files can contain reads from a whole sample (which can possibly
come from many different reference sequences) or can be demultiplexed
before alignment so that they contain reads from one reference sequence.


Options for FASTQ files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When giving FASTQ files on the command line, indicate the type of reads
(single-end, paired-end interleaved, or paired-end in separate files)
and the type of content (whole-sample or demultiplexed) via the option.
Each option has a "long form" (starting with ``--``) and a "short form"
(starting with ``-`` and comprising a single character), which are both
indicated in this table, separated by a ``/``:

========================== =============== =================
Endedness                  Whole Sample    Demultiplexed
========================== =============== =================
paired-end, separate files ``--fastqp/-x`` ``--dmfastqp/-X``
paired-end, interleaved    ``--fastqi/-y`` ``--dmfastqi/-Y``
single-end                 ``--fastqs/-z`` ``--dmfastqs/-Z``
========================== =============== =================

.. note::
    As with all options, the long and short forms are interchangable.
    For brevity, the rest of this manual uses only the short forms.

.. note::
    If using paired-end reads in separate FASTQ files, both files of 1st
    and 2nd mates are required. See :ref:`fastq_pair` for instructions.


Providing one FASTQ file (single-end or interleaved paired-end)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For example, to align a FASTQ of single-end reads from a whole sample,
use the option ``-z``::

    seismic align {refs.fa} -z {sample.fq.gz}

where ``{refs.fa}`` is the path to the FASTA file of reference sequences
and ``{sample.fq.gz}`` is the path to the FASTQ file of the sample.

For a FASTQ of paired-end, interleaved reads that were demultiplexed,
use the option ``-Y`` instead::

    seismic align {refs.fa} -Y {sample/ref.fq.gz}

where ``{sample/ref.fq.gz}`` is the path to the FASTQ file containing
reads from only one reference in the sample.


.. _fastq_pair:

Providing one pair of paired-end FASTQ files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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


Providing multiple FASTQ files or pairs or paired-end FASTQ files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

There are three ways to align multiple FASTQ files (or pairs) at once:

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


Options for Phred score encoding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Options for quality assessment with FastQC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, each FASTQ file is processed with `FastQC`_, both before and
after trimming, in order to find any potential problems.
FastQC can be disabled with the flag ``--no-fastqc``.
To enable automatic extraction of the zipped output files from FastQC,
add the flag ``--qc-extract``.


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
quality are trimmed with Cutadapt.
The default minimum quality is 25, which corresponds to a probability of
1 - 10 :sup:`-2.5` = 0.997 that the base call is correct.
(See :ref:`phred_encodings` for more details).
To change the quality threshold, use the option ``--min-phred``.

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

See :ref:`cli_align` for the full list of options that SEISMIC-RNA can
use with Cutadapt, and the `Cutadapt reference guide`_ for details on
each of these options.
These options should suffice for most users.
If you require a more customized adapter trimming workflow, you can trim
your FASTQ files outside of SEISMIC-RNA, then perform alignment within
SEISMIC-RNA using the option ``--no-cut`` to disable additional adapter
trimming.


Options for alignment with Bowtie 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bowtie 2 index files: automatic and pre-built
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA searches for a set of Bowtie 2 index files in the directory
where the FASTA file is.
A valid Bowtie 2 index comprises six files, all with the same name as
the FASTA file, with the extensions ``.1.bt2``, ``.2.bt2``, ``.3.bt2``,
``.4.bt2``, ``.rev.1.bt2``, and ``.rev.2.bt2``.
If this set of files exists, then SEISMIC-RNA uses it as the index.
Otherwise, it calls ``bowtie2-build`` to build an index in a temporary
directory, which is deleted after alignment finishes.
Indexing a small FASTA file takes several seconds; thus, it is generally
more convenient to build a temporary .
However, indexing a very large FASTA (e.g. a whole transcriptome) can
take hours, so it is advantageous to pre-build your index in the same
directory as the FASTA file, to save time in case you need to align more
than once.
You can pre-build an index with this command::

    bowtie2-build refs.fa {refs}

replacing ``{refs}`` with the path to and name of your FASTA file.
See the `Bowtie 2 Indexer manual`_ for more details.

.. note::
    If SEISMIC-RNA finds a pre-built Bowtie 2 index, then it does *not*
    verify that the index was actually built from the FASTA file of the
    same name.
    You must verify this yourself if using a pre-built index.
    You can assume the index is correct if you build it using the above
    command and avoid modifying or replacing the FASTA and index files.

Alignment modes: local and end-to-end
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

Filtering alignments
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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

Examining unaligned reads/pairs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For some reasons, including troubleshooting low alignment rates, it can
be helpful to examine the reads that did not align. Depending on the

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

Additional options for Bowtie 2
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

See :ref:`cli_align` for the full list of options that SEISMIC-RNA can
use with Bowtie 2, and the `Bowtie 2 manual`_ for details on each of
these options.
These options should suffice for most users.
If you require a more customized alignment workflow, you can align your
your FASTQ files outside of SEISMIC-RNA, then pass the resulting XAM
files into SEISMIC-RNA at the step :ref:`step_relate`.

Troubleshooting alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _low_align_rate:

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
3.  Reads that failed to align with Bowtie2 are written to the file
    ``out/{sample}/align/unaligned.`` Open the SAM file,
    process several unaligned reads randomly, and use `BLAST`_ to discern
    their origins, which can help in deducing what went wrong.


.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Cutadapt: https://cutadapt.readthedocs.io/en/stable/
.. _Cutadapt reference guide: https://cutadapt.readthedocs.io/en/stable/reference.html
.. _Bowtie 2 Indexer manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
.. _description of alignment modes in Bowtie 2: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment
.. _alignment score: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar
.. _section of the Bowtie 2 manual on alignment scores: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#valid-alignments-meet-or-exceed-the-minimum-score-threshold
.. _mapping quality: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mapping-quality-higher-more-unique
.. _view command in Samtools: https://www.htslib.org/doc/samtools-view.html
.. _Bowtie 2 manual for details on concordant/discordant alignments: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#concordant-pairs-match-pair-expectations-discordant-pairs-dont
.. _mixed mode: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mixed-mode-paired-where-possible-unpaired-otherwise
.. _overlap partially or completely, or dovetail: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other
.. _Bowtie 2 manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
