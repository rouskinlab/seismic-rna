********************************************************************************
seismic renumct
********************************************************************************


Purpose
================================================================================

``seismic renumct`` renumbers the positions in CT (connectivity table) files.
Use it when a CT file from an external source uses a different numbering scheme
than your reference sequence (e.g. starts at a position other than 1), and you
need the positions to match before using the file with SEISMIC-RNA.


Inputs
================================================================================

CT files or directories containing CT files
    Passed via ``--ct-pos-5 FILE FIRST`` (see Options below).
    See :doc:`/formats/data/ct`.


Outputs
================================================================================

Renumbered CT files.
By default, written to the output directory (``--out-dir``).
With ``--inplace``, the original files are overwritten.


Quick example
================================================================================

Renumber ``struct.ct`` so that position 1 in the file corresponds to position
34 in the reference (i.e. set the 5' position to 34)::

    seismic renumct --ct-pos-5 struct.ct 34


Options
================================================================================

``--ct-pos-5 FILE N``
    CT file (or directory of CT files) and the 1-based 5' position to assign.
    Repeat for multiple files.

``--inplace/--no-inplace``
    Overwrite the original files (default off).

    .. warning::

        ``--inplace`` is irreversible — back up your files first.

``--out-dir DIR`` (``-o``)
    Write renumbered files here (default ``./out``).
    Ignored when ``--inplace`` is set.

``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`/formats/data/ct` — CT file format
- :doc:`ct2db` — convert CT to dot-bracket format
- :doc:`db2ct` — convert dot-bracket to CT format
