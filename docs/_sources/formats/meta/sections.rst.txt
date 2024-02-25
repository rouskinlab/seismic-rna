
Metadata for Sections
------------------------------------------------------------------------

About metadata for sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for sections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Define and name sections of a reference sequence.
For references with many sections, it is more convenient, reproducible,
and traceable to define the sections in a file than on the command line,
using the option ``--sections-file`` (``-s``).
The sections file additionally permits giving each section a name.


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
  coordinate, 1-indexed and included in the section.
  An alternative way to specify the 5' end for amplicon-based samples is
  to type the forward primer in the field ``Forward Primer``.
  The primer must match the reference sequence exactly at exactly one
  location.
  The section's 5' end is placed downstream of the forward primer, with
  a gap whose length is set via the option ``--primer-gap``.
- The 3' end can be defined in the field ``3' End``, which is a numeric
  coordinate, 1-indexed and included in the section.
  An alternative way to specify the 3' end for amplicon-based samples is
  to type the reverse primer in the field ``Reverse Primer``.
  The reverse complement of the primer must match the reference sequence
  exactly at exactly one location.
  The section's 3' end is placed upstream of the reverse primer, with a
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

This text can be copied into a new CSV file.
