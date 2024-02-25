
List Positions Matching Criteria
--------------------------------------------------------------------------------

Background about lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use ``+listpos`` to list positions meeting specific criteria in tables.

Purpose of listing positions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

After you run ``seismic table``, you may want to find the positions that meet
certain criteria, such as the mutation rate being less/greater than a threshold.

How to list positions from a table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _listpos:

Command line for listing positions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Suppose you have already run ``seismic table``::

    seismic table out/sample/*/ref/sect

You can list the positions meeting certain criteria with ::

    seismic +listpos out/sample/table/ref/sect

This command will output one list file for each table file, with the same path
except ``table/`` is replaced by ``list/``.

You can filter the positions using options such as ``--max-fmut-pos 0.01``,
which will output all positions with a mutation rate of 0.01 or less.
Using ``--complement`` lists the complement of the positions; for example,
``--complement --max-fmut-pos 0.01`` would output all positions with a mutation
rate greater than 0.01.
