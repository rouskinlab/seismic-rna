
Relate: Compute relationships between references and aligned reads
--------------------------------------------------------------------------------

Relate: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _relate_refs:

Relate input file: Reference sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You need one file of reference sequences in FASTA format (for details on this
format, see :doc:`../../formats/data/fasta`).
If your file has characters or formatting incompatible with SEISMIC-RNA, then
you can fix it using the :doc:`../cleanfa` tool.

Relate input file: Alignment maps
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of alignment map files, each of which must be in SAM,
BAM, or CRAM (collectively, "XAM") format.
See :doc:`../../formats/data/xam` for more information.

.. note::
    The references in the FASTA file must match those to which the reads in the
    alignment map were aligned.
    Discrepancies can cause the Relate step to fail or produce erroneous output.
    You can assume that the references match if you use the same (unmodified)
    FASTA file for both the ``align`` and ``relate`` commands, or if you run
    both steps using the command ``seismic wf``.

Provide the alignment map files as a list after the FASTA file.
See :doc:`../inputs` for ways to list multiple files.
For example, to compute relation vectors for reads from ``sample-1`` aligned to
references ``ref-1`` and ``ref-2``, and from ``sample-2`` aligned to reference
``ref-1``, use the following command::

    seismic relate {refs.fa} sample-1/align/ref-1.cram sample-1/align/ref-2.cram sample-2/align/ref-1.cram

where ``{refs.fa}`` is the path to the file of reference sequences.

Relate: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Relate settings shared with alignment
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Because you can begin the SEISMIC-RNA workflow at ``seismic align`` or, if you
already have alignment map files, can begin at ``seismic relate``, these two
commands share the options ``--phred-enc``, ``--min-mapq``, ``--min-reads``, and
``--out-dir``.

Relate setting: Minimum Phred score
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In the Relate step, you can flag bases with low quality scores as ambiguous, as
if they were ``N``\s.
This step serves a purpose similar to that of quality trimming during the Align
step (see :ref:`quality_trimming`).
The difference is that quality trimming removes low-quality bases by shortening
reads from their ends, while the minimum quality score in the Relate step flags
low-quality bases located anywhere in the reads, while preserving read lengths.
See :ref:`relate_low_qual` for a more detailed description of how this works.

To set the minimum quality score, use ``--min-phred``.
The default is 25, meaning that base calls with a probabilities of at least
10\ :sup:`-2.5` = 0.3% of being incorrect are flagged as ambiguous.
(See :ref:`phred_encodings` for an explanation of quality scores.)
For example, if a ``T`` is called as a match with a quality score of 20, then it
would be flagged as possibly a match and possibly a subsitution to A, C, or G.

Relate setting: Ambiguous insertions and deletions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When insertions and deletions (indels) occur in repetitive regions, determining
which base(s) were inserted or deleted can be impossible due to the repetitive
reference sequence itself, even if the reads were perfectly free of errors.
To handle ambiguous indels, SEISMIC-RNA introduces a new algorithm that finds
all possible indels that could have produced the observed read (for details on
this algorithm, see :doc:`../algos/ambrel`).
This algorithm is enabled by default.
If you do not need to identify ambiguous indels, then you can disable this
algorithm with ``--no-ambrel``, which will speed up the Relate step at the cost
of reducing its accuracy on indels.

Relate setting: Batch size
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In the Relate step, you can divide up your data into batches to speed up the
analysis and reduce the amount of memory needed.
For an explanation of batching and how to use it, see :ref:`batches`.
You can specify batch size (in millions of base calls) using ``--batch-size``,
which is ``64.0`` (64 million base calls) by default.
Relate uses the batch size to calculate the number of reads in each batch.
The number of relationship bytes per batch, *B*, is the number of relationship
bytes per read, *L*, times the number of reads per batch, *N*:

*B* = *LN*

Since *L* is the length of the reference sequence and *B* is ``--batch-size``:

*N* = *B*/*L*

.. note::
    SEISMIC-RNA will aim to put exactly *N* reads in each batch but the last
    (the last batch can be smaller because it has just the leftover reads).
    If the reads are single-ended or were not aligned in `mixed mode`_, then
    every batch but the last will contain exactly *N* reads.
    If the reads are paired-ended and were aligned in `mixed mode`_, then
    batches may contain more than *N* reads, up to a maximum of 2\ *N* in the
    extreme case that only one read aligned in every mate pair.

Relate setting: Overhangs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When relating sequencing reads, SEISMIC-RNA will by default compute 
relationships for all base calls between the smallest 5' aligned position and 
the greatest 3' aligned position. 
This occurs independent of mate alignment orientation, and can result in the 
relating of base calls that fall outside the region between mate starts, for 
instance, if the 5' mate aligns to the negative strand at a position less than 
the 3' mate start position on the positive strand.
In certain rare circumstances, like when adapter trimming is inconsistent, or 
when using randomized adapter sequences during library preparation, this can 
result in SEISMIC-RNA calculating relationships for extraneous extensions.
The default behavior ``--overhangs`` can be disabled in favor of a more 
conservative approach ``--no-overhangs``, where only base calls greater than 
the 5' mate start and less than the 3' mate start positions 
(i.e. within the insert) are related.

Relate: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{sample}/relate/{ref}``, where
``{out}`` is the output directory, ``{sample}`` is the sample, and ``{ref}`` is
the name of the reference.

Relate output file: Batch of relation vectors
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Each batch of relation vectors contains a ``RelateBatchIO`` object and is saved
to the file ``relate-batch-{num}.brickle``, where ``{num}`` is the batch number.
See :doc:`../../data/relate/relate` for details on the data structure.
See :doc:`../../formats/data/brickle` for more information on brickle files.

Relate output file: Batch of read names
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Within each batch, the relate step assigns an index (a nonnegative integer) to
each read and writes a file mapping the indexes to the read names.
Each batch of read names contains a ``QnamesBatchIO`` object and is saved to the
file ``qnames-batch-{num}.brickle``, where ``{num}`` is the batch number.
See :doc:`../../data/relate/qnames` for details on the data structure.
See :doc:`../../formats/data/brickle` for more information on brickle files.

Relate output file: Reference sequence
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The relate step writes the reference sequence as a ``RefseqIO`` object to the
file ``refseq.brickle``.
See :doc:`../../data/relate/refseq` for details on the data structure.
See :doc:`../../formats/data/brickle` for more information on brickle files.

Relate output file: Relate report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA also writes a report file, ``relate-report.json``, that records the
settings you used for running the Relate step and summarizes the results.
See :doc:`../../formats/report/relate` for more information.

Relate: Troubleshoot and optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounted errors during the Relate step, then the most likely cause is
that the FASTA file or settings you used for the Relate step differ from those
that you used during alignment.

Insufficient reads in {file} ...
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This error means that you provided a SAM/BAM/CRAM file containing fewer reads
than the minimum number set by ``--min-reads`` (``-N``).
There are two common causes of this error:

- You ran ``seismic align`` and ``seismic relate`` separately (instead of with
  ``seismic wf``), and you used a larger value for ``--min-reads`` during the
  Relate step than the Align step.
  To check if this happened, open your report files from Align and Relate and
  see if the field "Minimum number of reads in an alignment map" has a larger
  value in the Relate report.
- You ran alignment outside of SEISMIC-RNA or obtained alignment map files from
  an external source, and some of the alignment maps have insufficient reads.

The solution for the problem is to ensure that you run ``seismic relate`` with
``--min-reads`` set to the minimum number of reads you actually want during the
Relate step.
As long as you do so, you may ignore error messages about insufficient reads,
since these messages just indicate that SEISMIC-RNA is skipping alignment maps
with insufficient reads, which is exactly what you want to happen.

Read {read} mapped with a quality score {score} ...
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This error means that a read inside an alignment file aligned with a mapping
quality lower than the minimum set by ``--min-mapq``.
There are two common causes of this error:

- You ran ``seismic align`` and ``seismic relate`` separately (instead of with
  ``seismic wf``), and you used a larger value for ``--min-mapq`` during the
  Relate step than the Align step.
  To check if this happened, open your report files from Align and Relate and
  see if the field "Minimum mapping quality to use an aligned read" has a larger
  value in the Relate report.
- You ran alignment outside of SEISMIC-RNA or obtained alignment map files from
  an external source, and some reads in the alignment maps have insufficient
  mapping quality.

The solution for the problem is to ensure that you run ``seismic relate`` with
``--min-mapq`` set to the minimum mapping quality you actually want during the
Relate step.
As long as you do so, you may ignore error messages about insufficient quality,
since these messages just indicate that SEISMIC-RNA is skipping reads with
with insufficient mapping quality, which is exactly what you want to happen.

Read {read} mapped to a reference named {name} ...
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This error means that a read inside an alignment file aligned to a reference
whose name does not match the name of the alignment file (minus the extension).
For example, if your alignment map file ``azure.cram`` contains a read that
aligned to a reference named ``cyan`` (instead of ``azure``), then you will get
this error message.

If you aligned the reads using ``seismic align`` or ``seismic wf``, then this
error should never occur (unless you renamed or modified the output files).
Otherwise, you can solve the problem by ensuring that

- Each alignment map file contains reads that aligned to only one reference.
- Each alignment map file is named (up to the file extension) the same as the
  one reference to which all of the reads aligned.

Relate crashes or hangs while producing few or no batch files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Most likely, your system has run out of memory.
You can confirm using a program that monitors memory usage (such as ``top`` in a
Linux/macOS terminal, Activity Monitor on macOS, or Task Manager on Windows).
If so, then rerun Relate with adjustments to one or both settings:

- Use smaller batches (with ``--batch-size``) to limit the size of each batch,
  at the cost of having more files with a larger total size.
- Use fewer processors (with ``--max-procs``) to limit the memory usage, at the
  cost of slower processing.

.. _mixed mode: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mixed-mode-paired-where-possible-unpaired-otherwise
