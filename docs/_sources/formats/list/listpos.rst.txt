
List of Positions
--------------------------------------------------------------------------------

List of Positions: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

========= ==== ========================================================= ========
Name      Type Description                                               Required
========= ==== ========================================================= ========
Reference str  Name of the reference sequence                            yes
Position  int  Position in the reference sequence                        yes
...       n/a  Additional field (ignored)                                no
========= ==== ========================================================= ========

List of Positions file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List of Positions as input file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``mask`` command accepts a List of Positions file via ``--exclude-file``.

List of Positions as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``+listpos`` command outputs a List of Positions file for each table file.

Example List of Positions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List of Positions as a pretty table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

========= ========
Reference Position
========= ========
refA            11
refA            16
refA            20
refB            11
refB            46
refB            83
========= ========

List of Positions as plain text
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
::

    Reference,Position
    refA,11
    refA,16
    refA,20
    refB,11
    refB,46
    refB,83

This text can be copied into a new CSV file.
