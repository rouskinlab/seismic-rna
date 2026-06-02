********************************************************************************
seismic ensembles
********************************************************************************


Purpose
================================================================================

``seismic ensembles`` automatically identifies independently folding segments
along a full-length RNA by finding positions whose mutation rates are correlated
with each other (correlated pairs) and then running ``seismic filter`` and
``seismic cluster`` on those segments.
Use it when you want an unbiased, data-driven discovery of structural
heterogeneity along a long RNA without defining regions manually.


Inputs
================================================================================

IDmut output directories or report files
    One or more IDmut output paths.
    See :doc:`/use/inputs`.


Outputs
================================================================================

All outputs go into subdirectories under the sample's output directory.
The step runs filter and cluster internally, producing the same output
structure as those steps for each discovered segment.


Quick example
================================================================================

Run ensembles discovery on a full-length RNA::

    seismic ensembles out/sample-1/idmut/long-rna


Options
================================================================================

Tiling
    ``--tile-length N`` (``-L``)
        Length of each tile in nucleotides (default 0 = 2× median read length).
    ``--tile-min-overlap F`` (``-O``)
        Minimum fractional overlap between adjacent tiles (default 0.5).
    ``--erase-tiles/--keep-tiles``
        Delete intermediate filter files from the tiling step (default: erase).

Correlated-pair detection
    ``--pair-fdr F``
        False discovery rate for detecting correlated pairs (default 0.05).
    ``--min-pairs N``
        Only cluster segments with at least N correlated pairs (default 2).

Segment length filters
    ``--min-cluster-length N``
        Only cluster segments with at least N positions (default 20).
    ``--max-cluster-length N``
        Only cluster segments with at most N positions (default 1200).

Gap handling
    ``--gap-mode {omit|insert|expand}``
        How to handle gaps between segments (default ``omit``).
    ``--threshold-divisor F``
        Divide the clustering threshold by F; increase for more (but less
        specific) modules (default 1.0).

Filter and Cluster options
    All ``seismic filter`` and ``seismic cluster`` options are accepted and
    forwarded to the internal runs.
    See :doc:`filter` and :doc:`cluster` for details.

Branches
    ``--branch NAME`` (``-b``)
        Write outputs under a branch name.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`filter` — run on each segment internally
- :doc:`cluster` — run on each segment internally
- :doc:`/use/branch`, :doc:`/use/parallel`
