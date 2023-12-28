
Dot-bracket (DB): RNA secondary structures
------------------------------------------------------------------------

DB file: Content format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a DB file, each RNA structure record comprises two parts:

- A header that contains the name of the structure.
- A body that contains the structure on one line (and for only the first
  structure, the RNA sequence on the line before).

See `DB format`_ for more information.

DB header line
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The first line of each structure record is the name of the structure.

DB body lines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For the first structure, the line after the header is the RNA sequence,
which can contain A, C, G, U, and N.

The structure is encoded on one line as a series of dots (``.``) and
brackets (``(``/``)``, ``[``/``]``, ``{``/``}``, ``<``/``>``).
A dot indicates that the base does not pair with another base.
A bracket indicates that the base pairs with the matching bracket.

DB file: Path format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DB file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA accepts the following extensions for dot-bracket files:

- ``.db`` (default)
- ``.dbn``
- ``.dot``

DB path parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Every output DB file is named after its mutational profile.

DB file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DB as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``fold`` command outputs a DB file for each predicted structure.

.. _DB format: https://rna.urmc.rochester.edu/Text/File_Formats.html
