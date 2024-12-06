
Join: Merge regions (horizontally) from the Mask or Cluster step
--------------------------------------------------------------------------------

Join: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Join input file: Mask/Cluster/Join report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of Mask or Cluster reports as inputs for the Join step.
You can also give Join report files, to join regions that were themselves made
by joining other regions.
See :doc:`../inputs` for ways to list multiple files.

.. note::
    SEISMIC-RNA will not double-count any of your original regions, even if
    they appear in multiple report files you are joining.
    It will just log a warning if it finds any regions given multiple times.

Mask, Cluster, and Join reports will be joined only if they share all of

- the top-level output directory, i.e. ``--out-dir`` (``-o``)
- the sample
- the reference
- whether they are clustered (i.e. regions from Mask will not be joined with
  regions from Cluster)

For each combination of these attributes, SEISMIC-RNA will produce a joined
region from all Mask/Cluster/Join reports with those attributes.
The original Mask/Cluster/Join report files will not be deleted or modified;
you will merely get a new Join report file for each joined region.

For example, if you ran the command ::

    seismic join -J {joined} {out}/{sample}/cluster/ref-1/reg-1 {out}/{sample}/cluster/ref-1/reg-2 {out}/{sample}/mask/ref-1/reg-1 {out}/{sample}/mask/ref-1/reg-2 {out}/{sample}/mask/ref-2/reg-1 {out}/{sample}/mask/ref-2/reg-2

where ``{out}`` is the path of your output directory and ``{joined}`` is the
name you want to give to each joined region, then you would get three new Join
reports representing the joined regions:

- ``{out}/{sample}/cluster/ref-1/{joined}/cluster-report.json``: made from
  ``{out}/{sample}/cluster/ref-1/reg-1/cluster-report.json`` and
  ``{out}/{sample}/cluster/ref-1/reg-2/cluster-report.json``
- ``{out}/{sample}/mask/ref-1/{joined}/mask-report.json``: made from
  ``{out}/{sample}/mask/ref-1/reg-1/mask-report.json`` and
  ``{out}/{sample}/mask/ref-1/reg-2/mask-report.json``
- ``{out}/{sample}/mask/ref-2/{joined}/mask-report.json``: made from
  ``{out}/{sample}/mask/ref-2/reg-1/mask-report.json`` and
  ``{out}/{sample}/mask/ref-2/reg-2/mask-report.json``

To join all valid combinations of Mask/Cluster/Join reports in ``{out}`` into
regions named ``{joined}``, you can use the command::

    seismic join -J {joined} {out}

Join: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Join setting: Name of joined region
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can choose the name of your joined region(s) using ``--joined`` (``-J``).
If you omit this option in ``seismic join``, then it will default to ``joined``.
If you omit this option in ``seismic wf``, then the Join step will not run.

Join setting: Clusters to join
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When joining cluster reports specifically, you may specify which cluster of each
region should be joined.
For example, if you are joining three regions, each of which you have clustered
up to order 2, then you may want cluster 1 of the joined region to comprise

- cluster 2 of region 1
- cluster 1 of region 2
- cluster 2 of region 3

and conversely, cluster 2 of the joined region to comprise

- cluster 1 of region 1
- cluster 2 of region 2
- cluster 1 of region 3

You can specify the clusters to join using ``--join-clusts`` (``-j``) followed
by a file of the clusters to join.
See :doc:`../../formats/meta/joined` for more information on this file.

Join: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directories ``{out}/{sample}/mask/{ref}/{joined}``
and/or ``{out}/{sample}/cluster/{ref}/{joined}``, where ``{out}`` is the output
directory, ``{sample}`` is the sample name, ``{ref}`` is the reference name,
and ``{joined}`` is the joined region name.

Join output file: Join report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA writes a Join report file called either ``mask-report.json`` or
``cluster-report.json``, depending on whether the joined regions came from Mask
or Cluster reports.
See :doc:`../../formats/report/join` for more information.
The file is named ``mask-report.json`` or ``cluster-report.json`` not because
its contents are identical to those of a Mask/Cluster report file (they aren't)
but because SEISMIC-RNA can more easily use them interchangably when they have
the same file names.
You can give Join reports to the Table step just as you would Mask or Cluster
reports.

Join: Troubleshoot and optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Joined region ... got duplicate regions ...
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This warning means that an original (unjoined) region appeared more than once
in the report files you are joining.

For example, suppose that you join ``region-A`` and ``region-B``::

    seismic join -J region-AB out/sample/mask/ref/region-A out/sample/mask/ref/region-B

Then you try to join ``region-A`` with the joined region ``region-AB``::

    seismic join -J region-AAB out/sample/mask/ref/region-A out/sample/mask/ref/region-AB

This second command will warn that ``region-A`` is duplicated because it
appears in both the report files for ``sample-A`` and ``pool-1``.

If you get this warning, then you should check your Pool report file to ensure
it contains all the samples you want and none that you don't.

Overwriting ... would cause data loss
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This error means that you attempted to create a joined region with the same
name as an existing non-joined region while using ``--force``, e.g. ::

    seismic join --force -J region-A out

if ``region-A`` already exists.

Doing so would overwrite the Mask/Cluster report for the original, non-joined
region, making it unusable.
To prevent data loss, the Join step refuses to overwrite Mask/Cluster reports,
even with ``--force``.
