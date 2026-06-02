********************************************************************************
seismic migrate
********************************************************************************


Purpose
================================================================================

``seismic migrate`` updates a SEISMIC-RNA output directory produced by an older
version so that it works with the current version.
The current migration handles version 0.24 → 0.25, which renamed two steps:

- ``relate/`` directories and files → ``idmut/``
- ``mask/`` directories and files → ``filter/``

and updated the JSON fields in report files to match renamed options.

The input directory is **never modified**.
``seismic migrate`` copies the entire output tree and updates the copy.


Inputs
================================================================================

Output directory from SEISMIC-RNA 0.24
    Pass the top-level output directory (the one whose immediate children
    are per-sample directories).
    Exactly one path may be supplied.


Outputs
================================================================================

A full copy of the input tree, with all renames and JSON field updates applied.
Written to ``--out-dir`` (default ``./out``).
The output directory must not exist before running the command.


Quick example
================================================================================

Migrate the output directory ``out/`` into a new directory ``out-0.25/``::

    seismic migrate out -o out-0.25

After the command succeeds, use ``out-0.25`` in place of ``out``.


Options
================================================================================

``--out-dir DIR`` (``-o``)
    Path for the migrated output (must not exist; default ``./out``).
``--num-cpus N``
    Parallel workers; one per sample directory (default: all CPUs).

The auto-generated :doc:`/cli` lists every option with its current default.


Caveats
================================================================================

- The output directory must not exist before the command runs.
- The output path must differ from the input path.
- Only the 0.24 → 0.25 migration is supported.
  Outputs from 0.23 or earlier may not migrate correctly.
- If ``idmut/{ref}/refseq.brickle`` is missing for a reference, all masked
  positions are attributed to G (the U list is left empty).


Common errors
================================================================================

*seismic migrate can process 1 directory at a time, but got N*
    You supplied more than one positional argument.

*<dir> contains no sample directories*
    You passed a sample directory or project root instead of the top-level
    output directory.

*the output directory must not exist*
    Delete or rename the existing directory, or choose a different path.


See also
================================================================================

- :doc:`/use/workflow/idmut` — the renamed IDmut step (formerly Relate)
- :doc:`/use/workflow/filter` — the renamed Filter step (formerly Mask)
