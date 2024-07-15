
Metadata for Joined Clusters
------------------------------------------------------------------------

About metadata for joined clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for joined clusters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Define which clusters for each order to join with ``seismic join``.
For example, suppose that you have clustered two sections separately,
each up to order 3, and now you want to join those clustered sections.
You could simply have cluster *i* of the joined section contain cluster
*i* of each individual section.
However, since the number of each cluster is arbitrary, you may want to
join clusters in other ways.
For example, you may want to join the clusters at order 3 like this:

- cluster 1 of the joined section contains

  - cluster 1 of section 1
  - cluster 3 of section 2

- cluster 2 of the joined section contains

  - cluster 3 of section 1
  - cluster 2 of section 2

- cluster 3 of the joined section contains

  - cluster 2 of section 1
  - cluster 1 of section 2

You can specify the cluster from each section that should go into each
joined cluster in the joined clusters metadata file.

.. note::
    When joining *n* sections, each clustered at order *k*, the number
    of ways to join those sections is *k*\ :sup:`n`.


Fields of metadata for joined clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======= ==== ============================================================== ========
Name    Type Description                                                    Required
======= ==== ============================================================== ========
Order   int  Order of the joined section                                    yes
Cluster int  Cluster of the joined section                                  yes
...     int  Cluster of the individual section to put in the joined section no
======= ==== ============================================================== ========

Notes about metadata for joined clusters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Fields are case-sensitive and must be in the file's first line.
- In the columns "Order" and "Cluster", the orders must start at 1 and
  increase with no gaps; and for each order, every cluster number from 1
  to the order must be given in increasing order on its own row.
- Every column after the first two ("Order" and "Cluster") should be
  named after one of the individual sections.
- In each section column, the number indicates which cluster from that
  section will be joined into the joined cluster whose order and number
  are given in the "Order" and "Cluster" columns.
  The number must be an integer between 1 and the order; and for each
  order, the cluster numbers may not be repeated.

Example metadata for joined clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Metadata for joined clusters as a pretty table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

===== ======= ========= ===========
Order Cluster mysection yoursection
===== ======= ========= ===========
    1       1         1           1
    2       1         1           2
    2       2         2           1
    3       1         1           3
    3       2         3           2
    3       3         2           1
===== ======= ========= ===========

Metadata for joined clusters as plain text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
::

    Order,Cluster,mysection,yoursection
    1,1,1,1
    2,1,1,2
    2,2,2,1
    3,1,1,3
    3,2,3,2
    3,3,2,1

This text can be copied into a new CSV file.
