
Metadata for Joined Clusters
----------------------------

About metadata for joined clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for joined clusters
"""""""""""""""""""""""""""""""""""""""

Define which clusters to join with ``seismic join`` for each number of
clusters (K).
For example, suppose that you have clustered two regions separately,
each into up to 3 clusters (K = 3), and now you want to join those
clustered regions.
You could simply have cluster *i* of the joined region contain cluster
*i* of each individual region.
However, since the number of each cluster is arbitrary, you may want to
join clusters in other ways.
For example, at K = 3 you may want to join the clusters like this:

- cluster 1 of the joined region contains

  - cluster 1 of region 1
  - cluster 3 of region 2

- cluster 2 of the joined region contains

  - cluster 3 of region 1
  - cluster 2 of region 2

- cluster 3 of the joined region contains

  - cluster 2 of region 1
  - cluster 1 of region 2

You can specify the cluster from each region that should go into each
joined cluster in the joined clusters metadata file.

.. note::
    When joining *n* regions, each clustered into *K* clusters, the number
    of ways to join those regions is *K*\ :sup:`n`.


Fields of metadata for joined clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======= ==== ============================================================ ========
Name    Type Description                                                  Required
======= ==== ============================================================ ========
K       int  Number of clusters (K) of the joined region                  yes
Cluster int  Cluster of the joined region                                 yes
...     int  Cluster of the individual region to put in the joined region no
======= ==== ============================================================ ========

Notes about metadata for joined clusters
""""""""""""""""""""""""""""""""""""""""

- Fields are case-sensitive and must be in the file's first line.
- In the columns "K" and "Cluster", the values of K must start at 1 and
  increase with no gaps; and for each K, every cluster number from 1
  to K must be given in increasing order on its own row.
- Every column after the first two ("K" and "Cluster") should be
  named after one of the individual regions.
- In each region column, the number indicates which cluster from that
  region will be joined into the joined cluster whose K and number
  are given in the "K" and "Cluster" columns.
  The number must be an integer between 1 and K; and for each
  K, the cluster numbers may not be repeated.

Example metadata for joined clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Metadata for joined clusters as a pretty table
""""""""""""""""""""""""""""""""""""""""""""""

= ======= ======== ==========
K Cluster myregion yourregion
= ======= ======== ==========
1 1       1        1
2 1       1        2
2 2       2        1
3 1       1        3
3 2       3        2
3 3       2        1
= ======= ======== ==========

Metadata for joined clusters as plain text
""""""""""""""""""""""""""""""""""""""""""
::

    K,Cluster,myregion,yourregion
    1,1,1,1
    2,1,1,2
    2,2,2,1
    3,1,1,3
    3,2,3,2
    3,3,2,1

This text can be copied into a new CSV file.
