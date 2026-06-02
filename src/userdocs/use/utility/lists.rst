********************************************************************************
seismic list
********************************************************************************


Purpose
================================================================================

``seismic list`` scans per-position mutation rate tables and writes a list of
positions that pass specified quality criteria.
Use it to identify high-quality or low-quality positions, or to create an
exclusion list for ``seismic filter``'s ``--mask-pos-file`` option.


Inputs
================================================================================

Table CSV files or table output directories
    Per-position mutation rate tables from ``seismic table``.
    See :doc:`/use/inputs`.


Outputs
================================================================================

One list file per input table, written to ``{out}/{sample}/list/{ref}/{reg}/``.

``{step}-position-list.csv``
    One row per position that passed the filters.
    See :doc:`/formats/list/listpos`.


Quick example
================================================================================

List positions with at least 1000 informative reads::

    seismic list out/sample-1/filter/ref-1


Options
================================================================================

Position quality thresholds
    ``--min-ninfo-pos N``
        Include only positions with at least N informative reads (default 1000).
    ``--max-fmut-pos F``
        Include only positions with at most fraction F mutated reads (default 1.0).

Branches
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/list_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`/use/workflow/table` — produces the tables this step scans
- :doc:`/use/workflow/filter` — accepts the output as ``--mask-pos-file``
- :doc:`/formats/list/listpos`
- :doc:`/use/branch`, :doc:`/use/parallel`
