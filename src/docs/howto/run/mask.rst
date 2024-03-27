
Mask: Define mutations and sections to filter reads and positions
--------------------------------------------------------------------------------

Mask: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mask input file: Relate report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of Relate report files as inputs for the Mask step.
See :doc:`../inputs` for ways to list multiple files.

For example, to mask relation vectors of reads from ``sample-1`` related to
references ``ref-1`` and ``ref-2``, and from ``sample-2`` related to reference
``ref-1``, use the command ::

    seismic mask {out}/sample-1/relate/ref-1 {out}/sample-1/relate/ref-2 {out}/sample-2/relate/ref-1

where ``{out}`` is the path of your output directory from the Relate step.

To mask all relation vectors in ``{out}``, you can use the command ::

    seismic mask {out}

Mask: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mask setting: Define sections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can mask the full reference sequences or select specific sections.
The latter is useful for investigating small elements of longer sequences, such
as a 350 nt `IRES`_ within a 9,600 nt viral genome.
See :doc:`../sections` for ways to define sections.

Mask setting: Define mutations
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Mask step takes in relation vectors -- which encode relationships including
ambiguous mutations -- and outputs bit vectors, wherein each position in each
read has a binary, mutated/matched status.
For more information on relation vectors, see :doc:`../../data/relate/codes`.

Producing bit vectors requires deciding which types of relationships count as
mutations, which count as matches, and which count as neither.
You can choose which types of relationships to count as matches and mutations.
The default is to count all 4 types of matches (A→A, C→C, G→G, T→T) as matches
and all 12 types of substitutions (A→C, A→G, A→T, C→A, C→G, C→T, G→A, G→C, G→T,
T→A, T→C, T→G) as mutations, but not to count deletions and insertions (indels).
To count deletions and insertions as mutations, add ``--count-del`` and
``--count-ins``, respectively.

You can also choose to not count individual types of relationships, such as
substitutions from A to G (but still count every other type of substitution).
To ignore one type of relationship, add ``--discount-mut`` followed by a code
of two lowercase letters:

- The first letter is the base in the reference (``a``/``c``/``g``/``t``)
- The second letter is the base in the read (for substitutions) or ``d``/``i``
  (for deletions and insertions, respectively).

For example, to count all substitutions except A→G and all deletions except
of C, use ``--count-del --discount-mut ag --discount-mut cd``.

.. _mask_exclude:

Mask setting: Exclude positions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The first substep of masking is excluding pre-specified positions.
You can specify three types of positions to exclude.

Exclude positions with G and U bases
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

DMS methylates G and U much less than A and C residues under physiological
conditions [`Zubradt et al. (2017)`_], so positions with G or U bases are
generally excluded when DMS is the chemical probe.
Use ``--exclude-gu`` (default) and ``--include-gu`` to choose whether to use
G and U bases.

Exclude positions with poly(A) sequences
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Although DMS and SHAPE reagents do modify A residues that are not immobilized
by base pairing, stretches of consecutive A residues tend to have very low
mutation rates due to an artifact from the reverse transcriptases that are used
commonly for mutational profiling (including TGIRT-III and SuperScript II)
[`Kladwang et al. (2020)`_].
Thus, using poly(A) sequences in structural analyses can give erroneous results.
SEISMIC-RNA automatically excludes all positions within stretches of 5 or more
consecutive A residues.
You can customize this behavior with ``--exclude-polya`` followed by the minimum
length of poly(A) sequences to exclude.
To disable poly(A) exclusion, use ``--exclude-polya 0``.

Exclude arbitary positions
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can also exclude any arbitary positions from any reference sequence.
A common reason to exclude a position is if the base is modified endogenously
in a way that causes mutations during reverse transcription.
To exclude an arbitrary position, use ``--exclude-file`` followed by a file of
all positions to exclude, in :doc:`../../formats/list/listpos` format.

Mask setting: Filter reads
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The second substep of masking is filtering reads.
You can filter reads based on three criteria, in this order:

Filter reads by number of positions covering the section
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can require every read to contain a minimum number of bases in the section
(i.e. set a minimum coverage) using ``--min-ncov-read`` followed by the minimum
coverage.
The minimum coverage must be at least 1 because reads that do not cover the
section at all should always be filtered out.
Note that this filter considers only positions that were not pre-excluded (see
:ref:`mask_exclude`).

Filter reads by fraction of informative positions
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can set a limit on the minimum information in each read, as a fraction of
the number of non-excluded positions in the read, using ``--min-finfo-read``
followed by the minimum fraction of informative positions.
For example, to require 95% of the non-excluded positions in the read to be
informative, use ``--min-finfo-read 0.95``.
Note that the denominator of this fraction is the number of bases in the read
that cover the section; it is not just the length of the section or of the read.

Filter reads by fraction of mutated positions
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Rarely, a read may have an excessive number of mutations, possibly because it
underwent template switching during reverse transcription or misaligned during
the Align step.
You can set a limit to the fraction of mutated positions in the read using
``--max-fmut-read``.
For example, using the default limit of 10%, a read with 121 informative and
15 mutated positions would have a mutated fraction of 15 / 121 = 12% and be
discarded, but a read with 121 informative and 10 mutated positions would have
a mutated fraction of 8% and be kept.
Using ``--max-fmut-read 1.0`` disables filtering by fraction mutated.

Filter reads by space between mutations
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Reads with closely spaced mutations are very underrepresented in mutational
profiling data, presumably because reverse transcripases struggle to read
through closely spaced pairs of modifications [`Tomezsko et al. (2020)`_].
Therefore, the data are biased towards reads without closely spaced mutations,
which would skew the mutation rates.
However, SEISMIC-RNA can correct the bias: first by removing any reads that
did happen to have mutations close together, then calculating the mutation
rates without such reads, and inferring what the mutation rates would have
been if no reads had dropped out.

The correction for observer bias is most important for finding alternative
structures and (to minimize surprises) does not run by default.
You can correct observer bias using ``--min-mut-gap`` followed by the minimum
number of non-mutated bases that must separate two mutations; reads with any
pair of mutations closer than this gap are discarded.
If you correct for observer bias, then we recommend using ``--min-mut-gap 3``,
based on our previous findings in `Tomezsko et al. (2020)`_.

Mask setting: Filter positions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The third substep of masking is filtering positions.
You can filter positions based on two criteria, in this order:

Filter positions by number of informative reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Estimating the fraction of mutated reads at a given position requires a large
number of reads so that the uncertainty (i.e. error bars) is much smaller than
the fraction of mutated reads.
The default minimum number of informative reads is 1000, which we have found
to yield a reasonably small uncertainties in the mutation fraction.
You can specify the minimum number of informative reads at each position using
``--min-ninfo-pos``.
We discourage going below 1000 reads unless you have multiple replicates, the
total number of informative reads at the position among all replicates is at
least 1000, and the mutation rates of the replicates correlate with a Pearson
or Spearman coefficient of at least 0.95.

Filter positions by fraction of mutated reads
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Mutational profiling generally yields fractions of mutated reads up to 0.3.
Positions with fractions of mutated reads that exceed 0.5 are likely to be
mutated for some reason other than chemcial probing, such as misalignment
(especially when two or more reference sequences are very similar), an
endogenous RNA modification (if the RNA came from cells), a mistake in the
template DNA (if the RNA was transcribed *in vitro*), or a mistake in the
reference sequence.
Thus, SEISMIC-RNA discards positions with a fraction of mutated reads greater
than 0.5, by default.
You can set the maximum fraction of mutated reads using ``--max-fmut-pos {f}``.

Mask: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{sample}/mask/{ref}/{sect}``,
where ``{out}`` is the output directory, ``{sample}`` is the sample, ``{ref}``
is the reference, and ``{sect}`` is the section.

Mask output file: Batch of masked reads
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Each batch of masked reads contains a ``MaskBatchIO`` object and is saved to the
file ``mask-batch-{num}.brickle``, where ``{num}`` is the batch number.
See :doc:`../../data/mask/mask` for more information on the data structure.
See :doc:`../../formats/data/brickle` for more information on brickle files.

Mask output file: Mask report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA also writes a report file, ``mask-report.json``, that records the
settings you used for running the Mask step and summarizes the results, such as
which and how many positions and reads were filtered out for each reason.
See :doc:`../../formats/report/mask` for more information.

Mask: Troubleshoot and optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _mask_too_many_reads:

Too many reads are filtered out
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In the Mask report file, check the settings for filtering reads and the number
of reads removed by each filter.

- If the settings appear too strict, then rerun the Mask step using new settings
  that would keep more reads, such as a lower value for ``--min-finfo-read`` or
  ``--min-mut-gap`` or a higher value for ``--max-fmut-read``.
- If you are losing too many reads for having too few informative positions,
  then also double check the 5' and 3' ends of the section over which you are
  masking and ensure that the section is not too long compared to your reads.
- If you are losing too many reads for having too many mutations, or mutations
  that are too close together, then there may be a problem with the data quality
  that is causing excessive mutations, such as

  - Your RNA was low-quality, contained many endogenous modififications that
    caused mutations during RT, or did not have the sequence you expected.
  - Your sequencing run gave low-quality base calls (check the FastQC reports)
    that you did not trim (in Align) or flag as ambiguous (in Relate).
  - You aligned to reference sequences that differ from the actual RNA.
  - Many reads misaligned (possibly because your FASTA file has several similar
    sequences), and your mapping quality filter did not remove misaligned reads.
  - In the Mask step, you did not pre-exclude problematic positions, such as
    sites of endogenous RNA modifications.

Too many positions are filtered out
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

In the Mask report file, check the settings for filtering positions and the
number of positions removed by each filter.

- If the settings appear too strict, then rerun the Mask step using new settings
  that would keep more positions, such as a lower value for ``--min-ninfo-pos``
  or a higher value for ``--max-fmut-pos``.
- If you are losing too many positions for having too few informative reads,
  then there are three likely reasons:

  - Your sample was sequenced with insufficient depth or quality.
  - Your sample contained insufficient RNAs from this reference/section.
  - You lost too many reads during filtering; see :ref:`mask_too_many_reads`.

- If you are losing too many positions for having too many mutations, then there
  may be a problem with the data quality that is causing excessive mutations,
  such as

  - Your RNA was low-quality, contained many endogenous modififications that
    caused mutations during RT, or did not have the sequence you expected.
  - Your sequencing run gave low-quality base calls (check the FastQC reports)
    that you did not trim (in Align) or flag as ambiguous (in Relate).
  - You aligned to reference sequences that differ from the actual RNA.
  - Many reads misaligned (possibly because your FASTA file has several similar
    sequences), and your mapping quality filter did not remove misaligned reads.

Mask crashes or hangs while producing few or no batch files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Most likely, your system has run out of memory.
You can confirm using a program that monitors memory usage (such as ``top`` in a
Linux/macOS terminal, Activity Monitor on macOS, or Task Manager on Windows).
If so, then you can either

- Use fewer processors (with ``--max-procs``) to limit the memory usage, at the
  cost of slower processing.
- Rerun Relate with smaller batches (with ``--batch-size``) to limit the size of
  each batch, at the cost of having more files with a larger total size.

.. _IRES: https://en.wikipedia.org/wiki/Internal_ribosome_entry_site
.. _Zubradt et al. (2017): https://doi.org/10.1038/nmeth.4057
.. _Kladwang et al. (2020): https://doi.org/10.1021/acs.biochem.0c00020
.. _Tomezsko et al. (2020): https://doi.org/10.1038/s41586-020-2253-5
