********************************************************************************
seismic export
********************************************************************************


Purpose
================================================================================

``seismic export`` converts SEISMIC-RNA table data into JSON files formatted
for import into external tools.
It is off by default in ``seismic wf``; enable it with ``--export``.


Inputs
================================================================================

Table CSV files or output directories
    Per-position mutation rate tables from ``seismic table``.
    See :doc:`/use/inputs`.


Outputs
================================================================================

``{sample}__webapp.json``
    All of the sample's data in one JSON file, written to the top of the
    output directory (one file per sample).


Quick example
================================================================================

Export all table output for ``sample-1``::

    seismic export out/sample-1/filter/ref-1


Options
================================================================================

Metadata
    ``--samples-meta FILE`` (``-S``)
        CSV file of sample metadata to embed in exported results.
        See :doc:`/formats/meta/samples`.
    ``--refs-meta FILE`` (``-R``)
        CSV file of reference metadata to embed in exported results.
        See :doc:`/formats/meta/refs`.

Positions
    ``--all-pos/--unmasked-pos``
        Export all positions (on) or only unmasked positions (off).
        Default: on.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`table` — produces the CSV tables this step exports
- :doc:`/formats/meta/samples`, :doc:`/formats/meta/refs`
- :doc:`/use/inputs`, :doc:`/use/parallel`
