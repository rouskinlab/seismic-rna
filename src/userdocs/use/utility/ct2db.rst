********************************************************************************
seismic ct2db
********************************************************************************


Purpose
================================================================================

``seismic ct2db`` converts RNA secondary structure files from CT (connectivity
table) format to dot-bracket format.
Use it when a downstream tool requires dot-bracket notation.


Inputs
================================================================================

CT files
    One or more ``.ct`` files.
    See :doc:`/formats/data/ct`.


Outputs
================================================================================

Dot-bracket files (``.db``) written alongside the input CT files.
See :doc:`/formats/data/db`.


Quick example
================================================================================

Convert a CT file to dot-bracket::

    seismic ct2db struct.ct


Options
================================================================================

``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`db2ct` — the reverse conversion
- :doc:`/formats/data/ct`, :doc:`/formats/data/db`
- :doc:`renumct` — renumber positions in CT files
