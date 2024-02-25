
Metadata for References
------------------------------------------------------------------------

About metadata for references
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of metadata for references
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Define metadata about reference sequences, such as the coordinates at
which barcodes are located.

Fields of metadata for references
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

========= ==== ========================================================= ========
Name      Type Description                                               Required
========= ==== ========================================================= ========
Reference str  Name of the reference sequence                            yes
...       n/a  Additional metadata of the reference                      no
========= ==== ========================================================= ========

Notes about metadata for references
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Fields are case-sensitive and must be in the file's first line.
- Additional fields with arbitrary names can be given.
  Their data types will be inferred as strings, integers, or floats.

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

This text can be copied into a new CSV file.
