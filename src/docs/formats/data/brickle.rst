
Brickle: Compressed Python Objects
------------------------------------------------------------------------

Brickle file: Content format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Brickle is not technically a data format.
It is rather a file containing a `pickled Python object`_ that has been
compressed with `brotli`_.
The name "brickle" can be a

- portmanteau of "brotli" and "pickle"
- suggestion that the file is supposed to contain a Python object that
  has been compressed into a dense brick of data
- deliberate misspelling of a `neighborhood in Miami, FL`_

As a brickle file can contain any Python object, the structure of the
object's data itself is of more interest than the file format.
See :doc:`../../data/index` for information on the data structures.

Brickle file: Path format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Brickle file extensions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA accepts the following extensions for brickle files:

- ``.brickle``

Brickle path parsing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Brickle files use a variety of path structures throughout SEISMIC-RNA.
They are most commonly used to store batches of data, in which case the
file name follows the format ``{btype}-batch-{bnum}.brickle``, where
``{btype}`` is the type of the batch and ``{bnum}`` is the number of the
batch (the first batch in each set is numbered 0).

Brickle file: Uses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Brickle as input file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Brickle files are not explicitly passed as input files via the command
line, but when a report file that is associated with brickle files is
passed as an input file, the brickle files will be read.
This process happens during the following commands:

- ``all``
- ``mask``
- ``cluster``
- ``table``

Brickle as output file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- The ``relate``, ``mask``, and ``cluster`` commands all output brickle
  files containing data along with their report files.

Brickle as temporary file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Brickle files are not used as temporary files by any command.

.. _`pickled Python object`: https://docs.python.org/3/library/pickle.html
.. _brotli: https://brotli.org/
.. _neighborhood in Miami, FL: https://en.wikipedia.org/wiki/Brickell
