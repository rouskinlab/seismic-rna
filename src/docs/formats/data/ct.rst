
Connectivity Table (CT): RNA secondary structures
------------------------------------------------------------------------

CT file: Content format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a CT file, each RNA structure record comprises two parts:

- A header that contains the name and length of the structure.
- A body that contains the structure.

See `CT format`_ for more information.

CT header line
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Every header line must start with the length of the structure, delimited
by whitespace.
Thereafter, from the first non-whitespace character to the end of the
line is the title of the structure.

CT body lines
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The remaining lines in each record encode the structure.
SEISMIC-RNA can parse sequences that obey the following rules:

- Every body line comprises six whitespace-delimited fields:

  1. Index in the section
  2. Base
  3. Previous (index minus 1), or 0 if the first index
  4. Next (index plus 1), or 0 if the last index
  5. Index of the pairing partner, or 0 if unpaired
  6. Position in the reference sequence

- The indexes must be contiguous positive integers beginning at 1.
- The base must be A, C, G, U, or N.
- The positions must be contiguous positive integers.

CT file: Path format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CT file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA accepts the following extensions for CT files:

- ``.ct``

CT path parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Every output CT file is named after its mutational profile, but names of
input CT files are not used for any purpose.

CT file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CT as input file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Structures for these commands must be input as CT files:

- ``seismic graph roc``
- ``seismic +renumct``

CT as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``fold`` command outputs a CT file for each predicted structure.
- The ``+renumct`` command outputs a renumbered CT file for each input
  CT file.

CT as temporary file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``fold`` command writes a temporary CT file from ``Fold``, which
  it then renumbers.

.. _CT format: https://rna.urmc.rochester.edu/Text/File_Formats.html
