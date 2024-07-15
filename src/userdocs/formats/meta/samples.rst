
Metadata for Samples
------------------------------------------------------------------------

About metadata for samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Define metadata about samples, such as the cell line or buffer in which
the RNA was probed, the type and concentration of the chemical probe,
and the temperature during chemical probing.


Fields of metadata for samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

============= ===== ================================================================ ========
Name          Type  Description                                                      Required
============= ===== ================================================================ ========
Sample        str   Name of the sample                                               yes
...           n/a   Additional metadata of the sample                                no
============= ===== ================================================================ ========

Notes about metadata for samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Fields are case-sensitive and must be in the file's first line.
- Additional fields with arbitrary names can be given.
  Their data types will be inferred as strings, integers, or floats.


Example metadata for samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example metadata for samples as a pretty table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

====== ==== ========== ========= ========================== ===== ============= =========== ======== ===
Sample User Date       Condition System                     Probe Concentration Temperature Duration ASO
====== ==== ========== ========= ========================== ===== ============= =========== ======== ===
noaso  Hui  2023-09-19 in vitro  sodium cacodylate (300 mM) DMS            10.5         310      240   0
aso1   Hui  2023-09-19 in vitro  sodium cacodylate (300 mM) DMS            10.5         310      240   1
aso2   Hui  2023-09-19 in vitro  sodium cacodylate (300 mM) DMS            10.5         310      240   2
ut     Hui  2023-09-19 in vitro  sodium cacodylate (300 mM) DMS             0.0         310      240   0
====== ==== ========== ========= ========================== ===== ============= =========== ======== ===

Example metadata for samples as plain text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
::

    Sample,User,Date,Condition,System,Probe,Concentration,Temperature,Duration,ASO
    noaso,Hui,2023-09-19,in vitro,sodium cacodylate (300 mM),DMS,10.5,310,240,0
    aso1,Hui,2023-09-19,in vitro,sodium cacodylate (300 mM),DMS,10.5,310,240,1
    aso2,Hui,2023-09-19,in vitro,sodium cacodylate (300 mM),DMS,10.5,310,240,2
    ut,Hui,2023-09-19,in vitro,sodium cacodylate (300 mM),DMS,0,310,240,0

This text can be copied into a new CSV file.
