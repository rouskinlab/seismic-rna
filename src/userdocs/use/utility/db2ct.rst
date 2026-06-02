********************************************************************************
seismic db2ct
********************************************************************************


Purpose
================================================================================

``seismic db2ct`` converts RNA secondary structure files from dot-bracket
format to CT (connectivity table) format.
Use it when a downstream tool requires CT notation.


Inputs
================================================================================

Dot-bracket files
    One or more ``.db`` files.
    See :doc:`/formats/data/db`.


Outputs
================================================================================

CT files (``.ct``) written alongside the input dot-bracket files.
See :doc:`/formats/data/ct`.


Quick example
================================================================================

Convert a dot-bracket file to CT::

    seismic db2ct struct.db


Options
================================================================================

``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`ct2db` — the reverse conversion
- :doc:`/formats/data/ct`, :doc:`/formats/data/db`
- :doc:`renumct` — renumber positions in CT files
