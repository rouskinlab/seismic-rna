
Relate Batch
------------------------------------------------------------------------

Each batch of relation vectors is a ``RelateBatchIO`` object.

Relate batch: Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following attributes encode the relationships between each read and
each position in the reference sequence:

========= ============================================ ====================================================================
Attribute Data Type                                    Description
========= ============================================ ====================================================================
``end5s`` ``numpy.ndarray[int]``                       array of the first position of the most upstream mate in each read
``mid5s`` ``numpy.ndarray[int]``                       array of the first position of the most downstream mate in each read
``mid3s`` ``numpy.ndarray[int]``                       array of the last position of the most upstream mate in each read
``end3s`` ``numpy.ndarray[int]``                       array of the last position of the most downstream mate in each read
``muts``  ``dict[int, dict[int, numpy.ndarray[int]]]`` array of the reads with each type of mutation at each position
========= ============================================ ====================================================================

.. note::
    The positions of the first and last bases in the reference sequence
    are defined to be 1 and the length of the sequence, respectively.

.. _relate_read_nums:

Relate batch: Structure of read numbers
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Each read, or pair of paired-end reads, is labeled with a non-negative
integer: 0 for the first read in each batch, and incrementing by 1 for
each subsequent read.
Within one batch, all read numbers are unique.
However, two different batches can have reads that share numbers.

Relate batch: Structure of 5' and 3' end positions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- ``end5s``, ``mid5s``, ``mid3s``, and ``end3s`` are all 1-dimensional
  ``numpy.ndarray`` objects.
- For any relate batch, ``end5s``, ``mid5s``, ``mid3s``, and ``end3s``
  all have the same length (which may be any integer ≥ 0).
- A read with index ``i`` corresponds to the ``i``\ th values of
  ``end5s``, ``mid5s``, ``mid3s``, and ``end3s``; denoted (respectively)
  ``end5s[i]``, ``mid5s[i]``, ``mid3s[i]``, and ``end3s[i]``.
- For every read ``i``:

  - 1 ≤ ``end5s[i]`` ≤ ``end3s[i]`` ≤ length of reference sequence
  - If paired-end and there is a gap of ≥ 1 nt between mates 1 and 2:

    - ``end5s[i]`` ≤ ``mid3s[i]`` < ``mid5s[i]`` ≤ ``end3s[i]``

  - Otherwise:

    - ``end5s[i]`` = ``mid5s[i]`` ≤ ``mid3s[i]`` = ``end3s[i]``

Relate batch: Structure of mutations
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``muts`` is a ``dict`` wherein

- each key is a position in the reference sequence (``int``)
- each value is a ``dict`` wherein

  - each key is a type of mutation (``int``, see :doc:`./codes` for more
    information)
  - each value is an array of the numbers of the reads that have the
    given type of mutation at the given position (``numpy.ndarray``)

Relate batch: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For example, suppose that the reference sequence is ``TCAGAACC`` and a
batch contains five paired-end reads, numbered 0 to 4:

==== ==== ============
Read Mate Alignment
==== ==== ============
0    1    ``_CAG____``
0    2    ``_____AGC``
1    1    ``___GTA__``
1    2    ``TCT_____``
2    1    ``____AAC_``
2    2    ``_CA_____``
3    1    ``TAAGT___``
3    2    ``______CC``
4    1    ``__AGA___``
4    2    ``___GA-C_``
Ref       ``TCAGAACC``
==== ==== ============

The positions, reads, and relationships can be shown explicitly as a
matrix (see :doc:`./codes` for information on the relationship codes):

==== === === === === === === === ===
Read 1   2   3   4   5   6   7   8
==== === === === === === === === ===
0    255 1   1   1   255 1   64  1
1    1   1   128 1   128 1   255 255
2    255 1   1   255 1   1   1   255
3    1   16  1   1   128 255 1   1
4    255 255 1   1   3   3   1   255
==== === === === === === === === ===

In a relate batch, they would be encoded as follows:

- ``end5s``: ``[2, 1, 2, 1, 3]``
- ``mid5s``: ``[4, 1, 3, 5, 3]``
- ``mid3s``: ``[6, 6, 5, 7, 7]``
- ``end3s``: ``[8, 6, 7, 8, 7]``
- ``muts``::

    {1: {},
     2: {16: [3]},
     3: {128: [1]},
     4: {},
     5: {3: [4], 128: [1, 3]},
     6: {3: [4]},
     7: {64: [0]},
     8: {}}

  Note that the numbers are shown here for visual simplicity as ``list``
  objects, but would really be ``numpy.ndarray`` objects.
