
Vectorization 
=============

The vectorization module determines whether each base in each sequencing read matches the reference RNA sequence.
It outputs a matrix in which each column is a base in the reference sequence, each row is a read,
and each element is a byte indicating the relationship between the bases in the read and reference sequences.
This matrix can then be used to calculate mutation frequencies or to cluster the reads to find alternative structures.


.. include:: examples.rst

.. include:: reference.rst

.. include:: IO.rst


.. _vector_algos:
.. include:: key_algos.rst

.. include:: edge_cases.rst

**Contributors**: Matthew Allan
