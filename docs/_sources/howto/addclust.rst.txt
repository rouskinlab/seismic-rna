
Add Clusters to an Already-Clustered Dataset
--------------------------------------------------------------------------------

Background about adding clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of adding clusters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

After running the Cluster step, you may want to continue clustering using more
clusters.
You could accomplish this by simply rerunning the Cluster step using a larger
``--max-clusters`` (``-k``).
But since Cluster always begins with 1 cluster, you would need to repeat all the
clusters that you had already run before being able to add more clusters, which
would waste your computer resources.
The Add Clusters tool lets you keep your existing clustering results and merely
append more clusters, which is faster and more efficient.

How to add clusters to a Cluster dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Command line for adding clusters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Suppose you have already run clustering with a maximum of 2 clusters (``-k 2``)
on a dataset::

    seismic cluster -k 2 {out}/{sample}/mask/{ref}/{sect}

To add more clusters, use ``seismic +addclust`` with the new maximum number of
clusters (e.g. 4) and the Cluster report as the input file::

    seismic +addclust -k 4 {out}/{sample}/cluster/{ref}/{sect}

This command will resume clustering with one more than the maximum number of
existing clusters (in this case, the previous maximum was 2, so ``+addclust``
will begin at 3).
The maximum number of clusters follows the same rules as in ``cluster`` (see
:ref:`cluster_max` for more information).
You do not need to specify any other settings for clustering, such as the number
of runs or the threshold for convergence: ``+addclust`` automatically uses the
same settings that you used originally in ``cluster``.

You can give any number of Cluster report files as inputs for Add Clusters.
See :doc:`./inputs` for ways to list multiple files.

The only new files produced are ``mus`` and ``props`` files for each additional
clustering run.
The Cluster report file, the batch files, and the counts file are all updated
in-place, as if you had run ``cluster`` the first time with the new ``-k``.
Because updating in-place means that there is a risk of data loss if any error
occurs during writing, the original files are all first backed up in a temporary
directory, which you can specify with ``--temp-dir`` (``-t``).
In case of a fatal error while updating the files, the original files will all
be restored from the backup, as if you had never run ``+addclust``.
