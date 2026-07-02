********************************************************************************
seismic clusterscan
********************************************************************************


Purpose
================================================================================

``seismic clusterscan`` clusters the **domains** that :doc:`filterscan`
detected along an RNA.
It runs the cluster step on each domain separately, inferring how many
alternative structures (clusters) that domain folds into and the proportion
of reads in each.
Use it after :doc:`filterscan` to resolve the structural heterogeneity within
each data-driven domain.


Inputs
================================================================================

Filterscan output directory or report file
    One or more ``filterscan`` output directories or ``filterscan-report.json``
    files.
    See :doc:`/use/inputs` for ways to select multiple inputs at once.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/clusterscan/{ref}/{reg}/``.

``clusterscan-report.json``
    Summary of the run, including the best number of clusters for each domain.
    See :doc:`/formats/report/clusterscan`.

The cluster step is also run on each domain, producing the same output as
:doc:`cluster` for every domain.


Quick example
================================================================================

Cluster the domains that ``filterscan`` found::

    seismic clusterscan out/sample-1/filterscan/long-rna


Options
================================================================================

Cluster options
    All :doc:`cluster` options are accepted and applied to each domain, such as
    ``--min-clusters``/``--max-clusters`` to bound the number of clusters and
    the various run-quality filters.

Branches
    ``--branch X`` (``-b``)
        Create a new branch: output results in ``{out}/{sample}/clusterscan_{branch}``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`filterscan` — detects the domains this step clusters
- :doc:`cluster` — run on each domain internally
- :doc:`/formats/report/clusterscan`
- :doc:`/use/inputs`, :doc:`/use/branch`, :doc:`/use/parallel`
