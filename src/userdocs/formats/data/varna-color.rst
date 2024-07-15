
VARNA Color: Color codes for VARNA
------------------------------------------------------------------------

If you draw RNA structures using the software `VARNA`_, then you can use
VARNA color files to color the bases by their mutation rates.

VARNA color file: Content format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each line of a VARNA color file corresponds to one position in the RNA
sequence you folded.
For each position that has data (i.e. is within the section from which
the data came and was not masked), the line has the normalized mutation
rate, which is a decimal number between 0 and 1 (inclusive).
For every position without data, the line just has the value -1.

VARNA color file: Path format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VARNA color file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA accepts the following extensions for VARNA color files:

- ``.txt``

VARNA color path parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Every output VARNA color file is named after its mutational profile,
plus ``__varna-color.txt`` at the end.

VARNA color file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VARNA color as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The Fold step outputs a VARNA color file for each predicted structure.

.. _VARNA: https://varna.lisn.upsaclay.fr/
