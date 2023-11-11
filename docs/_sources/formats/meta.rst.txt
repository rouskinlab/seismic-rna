
Metadata Formats
========================================================================

Using metadata files, information can be given about samples, reference
sequences, and sections.


Metadata for Samples
------------------------------------------------------------------------

About metadata for samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Define metadata about samples, such as the cell line or buffer in which
the RNA was probed, the type and concentration of the chemical probe,
and the temperature during chemical probing.

Commands that use metadata for samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- ``seismic all``
- ``seismic export``

File format of metadata for samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
comma-separated values (``.csv``)


Fields of metadata for samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

============= ===== ================================================================ ========
Name          Type  Description                                                      Required
============= ===== ================================================================ ========
Sample        str   Name of the sample                                               yes
User          str   Name(s) of the experimentalist(s)                                no
Date          str   Date of the experiment (YYYY-MM-DD)                              no
Condition     str   Experimental condition: "in vitro", "in cellulo", or "in vivo"   no
System        str   Buffer if in vitro, cell line if in cellulo, organism if in vivo no
Probe         str   Name of the chemical probe (e.g. "DMS", "1M7")                   no
Concentration float Concentration of the chemical probe (millimolar)                 no
Temperature   float Temperature of chemical probing (Kelvin)                         no
Duration      float Incubation time of chemical probing (seconds)                    no
...           n/a   Additional metadata of the sample                                no
============= ===== ================================================================ ========

Notes about metadata for samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Fields are case-sensitive and must be in the file's first line.
- Additional fields with arbitrary names can be given. Their data types
  will be inferred as either strings or floats.


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


Metadata for References
------------------------------------------------------------------------

About metadata for references
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for references
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Define metadata about reference sequences, such as the coordinates at
which barcodes are located.

Commands that use metadata for references
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- ``seismic all``
- ``seismic export``

File format of metadata for references
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
comma-separated values (``.csv``)

Fields of metadata for references
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

========= ==== ========================================================= ========
Name      Type Description                                               Required
========= ==== ========================================================= ========
Reference str  Name of the reference sequence                            yes
Barcode5  int  Coordinate of the barcode's 5' end (1-indexed, inclusive) no
Barcode3  int  Coordinate of the barcode's 3' end (1-indexed, inclusive) no
...       n/a  Additional metadata of the reference                      no
========= ==== ========================================================= ========

Notes about metadata for references
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Fields are case-sensitive and must be in the file's first line.
- Additional fields with arbitrary names can be given. Their data types
  will be inferred as either strings or floats.


Example metadata for references
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Metadata for references as a pretty table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

========= ======== ========
Reference Barcode5 Barcode3
========= ======== ========
MyFavRNA        10       17
OtherRef         1        6
========= ======== ========

Metadata for references as plain text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
::

    Reference,Barcode5,Barcode3
    MyFavRNA,10,17
    OtherRef,1,6


Metadata for Sections
------------------------------------------------------------------------

About metadata for sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for sections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Define and name sections of a reference sequence. For references with
many sections, it may be more convenient, reproducible, and/or trackable
to define the sections in a file than on the command line (using the
option ``--coords`` or ``--primers``). The sections file additionally
permits giving each section a name.

Commands that use metadata for sections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- ``seismic all``
- ``seismic mask``
- ``seismic fold``

File format of metadata for sections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

comma-separated values (``.csv``)


Fields of metadata for sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

============== ==== ========================================================= ========
Name           Type Description                                               Required
============== ==== ========================================================= ========
Section        str  Name of the section being defined                         yes
Reference      str  Name of the reference of which the section is part        yes
5' End         int  Coordinate of the section's 5' end (1-indexed, inclusive) no
3' End         int  Coordinate of the section's 3' end (1-indexed, inclusive) no
Forward Primer str  Sequence of the forward primer for the section            no
Reverse Primer str  Sequence of the reverse primer for the section            no
============== ==== ========================================================= ========

Notes about metadata for sections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Fields are case-sensitive and must be in the file's first line.
- The 5' end can be defined in the field ``5' End``, which is a numeric
  coordinate, 1-indexed and included in the section. An alternative way
  to specify the 5' end for samples prepared as amplicons is to type the
  forward primer in the field ``Forward Primer``. The primer must match
  the reference sequence exactly at exactly one location. The section's
  5' end is placed downstream of the forward primer, with an intervening
  gap whose length is set via the option ``--primer-gap``.
- The 3' end can be defined in the field ``3' End``, which is a numeric
  coordinate, 1-indexed and included in the section. An alternative way
  to specify the 3' end for samples prepared as amplicons is to type the
  forward primer in the field ``Reverse Primer``. The primer must match
  the reference sequence exactly at exactly one location. The section's
  3' end is placed upstream of the forward primer, with an intervening
  gap whose length is set via the option ``--primer-gap``.


Example metadata for sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Metadata for sections as a pretty table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

======= ========= ====== ====== ============== ==============
Section Reference 5' End 3' End Forward Primer Reverse Primer
======= ========= ====== ====== ============== ==============
5utr    MyFavRNA       1    103
cds     MyFavRNA     104   2368
3utr    MyFavRNA    2369   2695
thing1  OtherRef                ACCCGTAACTATCG TACAGGTCCGCATG
======= ========= ====== ====== ============== ==============

Metadata for sections as plain text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
::

    Section,Reference,5' End,3' End,Forward Primer,Reverse Primer
    5utr,MyFavRNA,1,103,,
    cds,MyFavRNA,104,2368,,
    3utr,MyFavRNA,2369,2695,,
    thing1,OtherRef,,,ACCCGTAACTATCG,TACAGGTCCGCATG
