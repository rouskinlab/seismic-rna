
Table: Count mutations for each read and position
--------------------------------------------------------------------------------

Table: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Table input file: Report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of report files from the Relate, Mask, and Cluster steps
as inputs for the Table step.
See :doc:`../inputs` for ways to list multiple files.

To tabulate all results in ``{out}``, you can use the command ::

    seismic table {out}

Table: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Table setting: Suppress per-position and per-read tables
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA outputs all possible tables by default.
To suppress per-position and per-read tables, use ``--no-table-pos`` and
``--no-table-read``, respectively.

Table: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{sample}/table/{ref}/{sect}``,
where ``{out}`` is the output directory, ``{sample}`` is the sample, ``{ref}``
is the reference, and ``{sect}`` is the section.
Each output file from Relate, Mask, and Cluster reports are prefixed with
``relate-``, ``mask-``, and ``clust-``, respectively.

Table output file: Relationships per position
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA outputs the number of relationships at each position to a CSV file
named ``{step}-per-pos.csv``, where ``{step}`` is ``relate``/``mask``/``clust``.

Table output file: Relationships per read
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA outputs the number of relationships in each read to a CSV file named
``{step}-per-read.csv.gz``, where ``{step}`` is ``relate``/``mask``/``clust``.

Table output file: Reads per cluster
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For Cluster reports, SEISMIC-RNA outputs the number of reads in each cluster to
a CSV file named ``clust-freq.csv``.
