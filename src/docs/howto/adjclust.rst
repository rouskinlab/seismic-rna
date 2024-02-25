
Add/Delete Orders to/from an Already-Clustered Dataset
--------------------------------------------------------------------------------

Background about adding/deleting orders of clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use ``+addclust`` and ``+delclust`` to adjust the maximum number of
clusters in (i.e. maximum order of) a dataset that you have already clustered
with ``seismic cluster``.

Purpose of adding clusters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

After you run ``seismic cluster``, you may wish that you had used a higher limit
for the number of clusters (i.e. maximum order).
You could increase the maximum order simply by rerunning ``seismic cluster``
using ``--force`` and a larger value for ``--max-clusters`` (``-k``).
Since ``seismic cluster`` always begins with 1 cluster, you would need to repeat
all the orders that you had already run before being able to add more -- a waste
of your computer resources.
A more efficient command to increase the maximum order is ``seismic +addclust``,
which merely tests higher orders of clustering and, if successful, appends them
to your existing results.

Purpose of deleting clusters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

After you run ``seismic cluster``, you may realize that higher orders are not
useful for your analysis or have unwanted artifacts.
If you simply ignored them, they would still waste space on your file system:
higher orders result in larger files of cluster batches and cluster tables and
more numerous files during graphing.
You could decrease the maximum order simply by rerunning ``seismic cluster``
using ``--force`` and a smaller value for ``--max-clusters`` (``-k``).
Since your clustering results would already exist -- you would merely want to
delete some of them -- any approach that would require rerunning clustering
would not be the most efficient.
A more efficient command to decrease the maximum order is ``seismic +delclust``,
which merely deletes higher orders of clustering from your results.

How to add/delete orders to/from a Cluster dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _addclust:

Command line for adding orders
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Suppose you have already run clustering with a maximum order of 3 (``-k 3``)::

    seismic cluster -k 3 {out}/{sample}/mask/{ref}/{sect}

If this dataset stopped at order 1 or 2 (because the BIC failed to decrease),
then you cannot add any higher orders.
If this dataset stopped at order 3 (because that was the limit), then you can
add higher orders using ``seismic +addclust``, specifying a new maximum order
(e.g. 5) and the Cluster report(s) as input file(s)::

    seismic +addclust -k 5 {out}/{sample}/cluster/{ref}/{sect}

This command will start clustering with one more than the previous maximum order
(in this case, the previous maximum was 3, so ``+addclust`` will begin at 4).
The rules for stopping based on the BIC are the same as in ``seismic cluster``
(see :ref:`cluster_max` for more information).
You do not need to specify any other settings for clustering, such as the number
of runs or the threshold for convergence: ``+addclust`` automatically uses the
same settings that you used originally in ``cluster``.

You can give any number of Cluster report files as inputs for ``+addclust``.
See :doc:`./inputs` for ways to list multiple files.

The only new files produced are ``mus`` and ``props`` files for each additional
clustering run.
The Cluster report file, the batch files, and the counts file are all updated
in-place, as if you had run ``cluster`` the first time with the new ``-k``.

.. _delclust:

Command line for deleting orders
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Suppose you have already run clustering with a maximum order of 3 (``-k 3``)::

    seismic cluster -k 3 {out}/{sample}/mask/{ref}/{sect}

You can delete orders above a limit using ``seismic +delclust``, specifying a
new maximum order (e.g. 2) and the Cluster report(s) as input file(s)::

    seismic +delclust -k 2 {out}/{sample}/cluster/{ref}/{sect}

This command will simply delete the results for higher orders (in this case, 3)
from the Cluster report file, the batch files, and the counts file.
It does not delete the ``mus`` and ``props`` files, since these files are small,
not used by any subsequent step, and not constrained by the maximum order.

You can give any number of Cluster report files as inputs for ``+addclust``.
See :doc:`./inputs` for ways to list multiple files.

The Cluster report file, the batch files, and the counts file are all updated
in-place, as if you had run ``cluster`` the first time with the new ``-k``.

.. note::
    Be careful with ``+delclust`` because it will delete all orders greater than
    ``max-clusters`` (``-k``).
    Deleted orders cannot be restored directly; to regenerate them, you would
    need to rerun clustering with ``+addclust``.
