********************************************************************************
seismic table
********************************************************************************


Purpose
================================================================================

``seismic table`` summarizes per-position and per-read mutation data into CSV
files that can be opened in spreadsheet software or used for plotting.
Run it after Filter, Cluster, or Join to get the mutation rates and read counts
that power downstream graphs and structure prediction.


Inputs
================================================================================

IDmut, Filter, Cluster, or Join output directories or report files
    One or more paths.
    See :doc:`/use/inputs`.


Outputs
================================================================================

Output files go into the same directory as the input report they summarize.

Per-position tables
    ``{step}-position-table.csv``
        One row per reference position; columns include mutation rate,
        informative read count, and per-type mutation counts.

Per-read tables
    ``{step}-read-table.csv``
        One row per read; columns include the number of mutations.

Per-cluster tables
    ``cluster-abundance-table.csv``
        For Cluster data, one row per cluster with its number of reads.


Quick example
================================================================================

Tabulate filter output for ``sample-1``::

    seismic table out/sample-1/filter/ref-1

Tabulate everything in the output directory::

    seismic table out/


Options
================================================================================

Toggle which tables to write
    ``--idmut-pos-table/--no-idmut-pos-table``
        Write per-position table for IDmut data (default on).
    ``--idmut-read-table/--no-idmut-read-table``
        Write per-read table for IDmut data (default off).
    ``--filter-pos-table/--no-filter-pos-table``
        Write per-position table for Filter data (default on).
    ``--filter-read-table/--no-filter-read-table``
        Write per-read table for Filter data (default on).
    ``--cluster-pos-table/--no-cluster-pos-table``
        Write per-position table for Cluster data (default on).
    ``--cluster-abundance-table/--no-cluster-abundance-table``
        Write per-cluster read-count table for Cluster data (default on).

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`filter`, :doc:`cluster`, :doc:`join` — produce the data this step
  summarizes
- :doc:`graph` — plots the tables produced here
- :doc:`fold` — uses per-position mutation rates to predict structure
- :doc:`/use/inputs`, :doc:`/use/parallel`
