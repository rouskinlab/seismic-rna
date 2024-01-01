
Workflow: Run all steps
--------------------------------------------------------------------------------

For convenience, ``seismic wf`` runs all steps of the main workflow:

.. image::
    wf.png

.. note::
    When you run the workflow with ``seismic wf``, the steps Cluster, Fold, and
    Export do not run by default, but you can turn them on with options:

    - To cluster, use ``--max-clusters`` (``-k``) followed by the maximum number
      of clusters (in ``seismic wf``, ``-k`` defaults to 0 -- no clustering).
    - To fold, use ``--fold``.
    - To export, use ``--export``.

SEISMIC-RNA will process every input file from the point that it enters the
workflow, following the gray arrays, until the end of the entire workflow.
For example, for each file of read sequences that you input with ``-x``, ``-y``,
or ``-z``, it will

1.  run Align on the FASTQ files of reads to yield alignment maps
2.  run Relate on the alignment maps to yield Relate reports
3.  run Mask on the Relate reports to yield Mask reports
4.  (optionally) run Cluster on the Mask reports to make Cluster reports
5.  run Table on the Relate, Mask, and (optionally) Cluster reports to yield
    table files
6.  (optionally) run Fold on the Mask and (optionally) Cluster table files to
    yield predicted RNA secondary structure models
7.  (optionally) run Export on the Relate, Mask, and (optionally) Cluster table
    files and (optionally) structure models to yield sample results
8.  run Graph on the Relate, Mask, and (optionally) Cluster tables and
    (optionally) structure models to yield graph files

By contrast, for each Mask table file you input, SEISMIC-RNA will only

1.  (optionally) run Fold on the table file to yield predicted RNA secondary
    structure models
2.  (optionally) run Export on the table file and (optionally) structure models
    to yield sample results
3.  run Graph on the table file and (optionally) structure models to yield graph
    files

Workflow: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can give any input file to ``seismic wf`` if you could give it to any
individual command in the workflow.

The only required input file is of the reference sequences, which is used by the
Align and Relate steps; see :ref:`relate_refs` for more information.
This file must be the first positional argument.

You can list additional input files (e.g. alignment maps, report files, table
files) as positional arguments after the reference sequences, using any of the
ways to list input files (see :doc:`../inputs`), for example ::

    seismic wf {refs.fa} sample1/refA.bam {out}/*/relate/refB/relate-report.json {out}/*/table/

where ``{refs.fa}`` is the path of your FASTA file of reference sequences and
``{out}`` is your top-level output directory.

You can put keyword input files (``-x``/``-y``/``-z``/``-X``/``-Y``/``-Z``) in
any locations amid positional arguments, for example ::

    seismic wf {refs.fa} -x sample4/ -y sample5.fq.gz -z sample6.fq.gz sample1/refA.bam {out}/*/relate/refB/relate-report.json {out}/*/table/

Workflow: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use any option in ``seismic wf`` if you could use it in any individual
command in the workflow.
Additionally, ``--fold``/``--no-fold`` and ``--export``/``--no-export`` control
whether the Fold and Export steps, respectively, will run in ``seismic wf``
(both are "no" by default).
Options can be mixed in any locations amid positional arguments, for example ::

    seismic wf --force {refs.fa} -x sample4/ --min-finfo-read 0.95 --min-ninfo-pos 100000 -k 3 {out}/*/table/ --fold

Workflow: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``seismic wf`` command creates all output files that would have been created
by running every step of the workflow individually.
All output files go into the directory that you specify with ``--out-dir``,
which is ``./out`` by default.

Workflow: Troubleshoot and optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Refer to the advice for the individual step you are troubleshooting/optimizing.
