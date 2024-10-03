
Export: Export a file of each sample for the `seismic-graph`_ web app
--------------------------------------------------------------------------------

Export: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Export input file: Table files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of table files of masked or clustered reads as inputs:

- ``mask-per-pos.csv``
- ``mask-per-read.csv``
- ``clust-per-pos.csv``
- ``clust-per-read.csv``
- ``clust-freq.csv``

See :doc:`../inputs` for ways to list multiple files.
(SEISMIC-RNA will just ignore ``relate-per-pos.csv``/``relate-per-read.csv.gz``
files.)

To export all data for all samples in the output directory ``{out}``, you could
use the command ::

    seismic export {out}

To export all data for sample ``{sample}`` in the output directory ``{out}``,
you could use the command ::

    seismic export {out}/{sample}

To export only masked data for reference ``{ref}`` in the output directory
``{out}``, you could use the command ::

    seismic export {out}/*/{ref}/*/mask*

Export: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Export setting: Metadata
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can add arbitrary metadata for each sample using ``--samples-meta`` (``-S``)
followed by a file containing :doc:`../../formats/meta/samples`.
You can add arbitrary metadata for each reference using ``--refs-meta`` (``-R``)
followed by a file containing :doc:`../../formats/meta/refs`.

SEISMIC-RNA expects you to include metadata, so if you do not, then it will log
a warning (but will still export the data without metadata).

Export: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Export output file: Sample results
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The file of sample results goes into the top-level output directory ``{out}``.
See :doc:`../../formats/data/sample` for information on this file.

.. _seismic-graph: https://rouskinlab.github.io/seismic-graph/
