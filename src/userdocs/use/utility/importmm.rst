********************************************************************************
seismic importmm
********************************************************************************


Purpose
================================================================================

``seismic importmm`` converts RNA Framework Mutation Map (MM) files into the
IDmut batch format used by SEISMIC-RNA.
Use it if you generated mutation data with RNA Framework or DRACO and want to
continue analysis (filter, cluster, fold, etc.) in SEISMIC-RNA without
re-aligning.


Inputs
================================================================================

Mutation Map (MM) files
    One or more ``.mm`` files from RNA Framework or DRACO.
    See :doc:`/use/inputs`.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/idmut/{ref}/``, mirroring the output of
``seismic idmut``.

``idmut-batch-{num}.brickle``
    Per-read, per-position mutation calls.

``idmut-report.json``
    Summary of settings and results.


Quick example
================================================================================

Import MM files for a sample named ``my-sample``::

    seismic importmm --sample my-sample data/*.mm


Options
================================================================================

Required
    ``--sample NAME`` (``-s``)
        Name to give the imported sample.

Quality filters
    ``--min-reads N`` (``-N``)
        Skip MM files with fewer than N reads (default 1000).
    ``--batch-size N``
        Maximum reads per output batch (default 65536).

Branches
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/idmut_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`/use/workflow/idmut` — the step whose output format this produces
- :doc:`/use/workflow/filter` — next step after importing
- :doc:`/formats/data/brickle`
- :doc:`/use/branch`, :doc:`/use/parallel`
