
Read Names Batch
------------------------------------------------------------------------

Each batch of relation vectors is a ``QnamesBatchIO`` object.

Read names batch: Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The read names are in the attribute ``names``, a ``numpy.ndarray`` with
data type ``str``.
The name of the read with number ``i`` (see :ref:`relate_read_nums` for
more information) is at index ``i`` of the ``names`` array.

Read names batch: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose a XAM file contains five reads:

=========== ===========================================
Read Number Read Name
=========== ===========================================
0           ``VL00355:66:AACYHM7M5:1:2106:54095:30496``
1           ``VL00355:66:AACYHM7M5:1:2303:34554:55259``
2           ``VL00355:66:AACYHM7M5:1:2411:30558:38352``
3           ``VL00355:66:AACYHM7M5:1:1206:41825:52949``
4           ``VL00355:66:AACYHM7M5:1:2506:10903:12605``
=========== ===========================================

The relate step would write a batch with the ``names`` attribute ::

    ["VL00355:66:AACYHM7M5:1:2106:54095:30496",
     "VL00355:66:AACYHM7M5:1:2303:34554:55259",
     "VL00355:66:AACYHM7M5:1:2411:30558:38352",
     "VL00355:66:AACYHM7M5:1:1206:41825:52949",
     "VL00355:66:AACYHM7M5:1:2506:10903:12605"]

Note that the attribute is shown as a ``list`` for visual simplicity,
but would really be a ``numpy.ndarray``.
