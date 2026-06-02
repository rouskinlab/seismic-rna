********************************************************************************
seismic pool
********************************************************************************


Purpose
================================================================================

``seismic pool`` merges IDmut data from multiple samples into a single pooled
sample.
Use it when you sequenced the same RNA across multiple sequencing runs or
libraries and want to analyze all reads together.
The original samples are not modified; pool creates a new virtual sample that
points to all of them.


Inputs
================================================================================

IDmut output directories or report files
    One or more paths pointing to IDmut (or previously pooled) output.
    See :doc:`/use/inputs`.

Pooled sample name
    First positional argument: the name you want to give the new pooled sample.


Outputs
================================================================================

All outputs go into ``{out}/{pooled-sample}/idmut/{ref}/``.

``idmut-report.json``
    Pool report listing the names of all contributing samples.
    See :doc:`/formats/report/pool`.
    This file is named ``idmut-report.json`` so that downstream steps
    (Filter, Table) accept it interchangeably with IDmut reports.


Quick example
================================================================================

Pool ``sample-1`` and ``sample-2`` into a sample named ``pooled``::

    seismic pool pooled out/sample-1/idmut out/sample-2/idmut

To pool every IDmut output in the output directory::

    seismic pool pooled out/


Options
================================================================================

Quality filters
    ``--min-pearson F``
        Skip pooling a pair of samples if their Pearson correlation is below F
        (default 0.0, i.e. no filter).
    ``--max-marcd F``
        Skip pooling a pair of samples if their mean arcsine distance exceeds F
        (default 1.0, i.e. no filter).

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

Warning: duplicate samples
    An original sample appears more than once among the inputs (e.g. because
    you included both a sample and a pool that already contains it).
    Check the resulting pool report to confirm it contains what you expect.

Error: would overwrite existing non-pooled sample
    You used a pooled sample name that already exists as a real sample.
    Choose a different name or remove the conflicting directory first.


See also
================================================================================

- :doc:`idmut` — produces the data this step pools
- :doc:`filter` — accepts pool output alongside idmut output
- :doc:`/formats/report/pool`
- :doc:`/use/inputs`, :doc:`/use/parallel`
