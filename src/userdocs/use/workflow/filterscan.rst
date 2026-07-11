********************************************************************************
seismic filterscan
********************************************************************************


Purpose
================================================================================

``seismic filterscan`` searches a full-length RNA for **domains**: regions that
appear to fold independently, revealed by positions whose mutations are
correlated with each other (correlated pairs).
It slides overlapping tiles along the RNA, runs the filter step on each tile,
finds pairs of positions that mutate together more often than expected by
chance, and groups those pairs into domains.
It then runs the filter step once more on each detected domain, so that the
domains are ready to be clustered.
Use it when you want an unbiased, data-driven way to discover structural
domains along a long RNA without defining regions by hand.
Follow it with :doc:`clusterscan` to cluster the domains it detects.


Inputs
================================================================================

IDmut output directory or report file
    One or more IDmut output directories or ``idmut-report.json`` files.
    See :doc:`/use/inputs` for ways to select multiple inputs at once.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/filterscan/{ref}/{reg}/``.

``filterscan-report.json``
    Summary of settings and results, including the coordinates of every
    detected domain.
    See :doc:`/formats/report/filterscan`.

``pairs.csv``
    The correlated pairs of positions found in the RNA.

``domains.csv``
    The coordinates (5' and 3' ends) of the detected domains.

``pairs_and_domains.html``
    An interactive plot of the correlated pairs and the domains built from them.

``confusion-matrix.csv``
    The per-pair counts and statistics used to decide which pairs are correlated.

The filter step is also run on each detected domain, producing the same output
as :doc:`filter` for every domain.


Quick example
================================================================================

Scan a full-length RNA for domains::

    seismic filterscan out/sample-1/idmut/long-rna


Options
================================================================================

Tiling
    ``--tile-length N`` (``-L``)
        Length of each tile in nucleotides (default 0 = 2× the median read length).
    ``--tile-min-overlap F`` (``-O``)
        Minimum fractional overlap between adjacent tiles (default 0.5).
    ``--erase-tiles/--keep-tiles``
        Delete the intermediate filter files from the tiling step (default: erase).

Correlated-pair detection
    ``--pair-fdr F``
        False discovery rate for calling a pair of positions correlated
        (default 0.05).
    ``--pair-distance-percentile F``
        Drop a correlated pair if its nearest surviving neighbor is farther than
        this percentile of all nearest-neighbor distances (default 95.0).
        Isolated pairs are treated as noise.
    ``--min-nearby-pairs N``
        Minimum number of other surviving pairs that must lie within the
        distance threshold for a pair to be kept (default 2).
        Values above 1 filter out small coincidental clusters of noise pairs.

Domain length filters
    ``--min-pairs N``
        Keep only domains with at least this many correlated pairs (default 2).
    ``--min-cluster-length N``
        Keep only domains with at least this many positions (default 20).

Gap handling
    ``--gap-mode {omit|insert|expand}``
        What to do with the gaps between domains (default ``omit``):
        ``omit`` leaves the gaps out, ``insert`` turns each gap into its own
        domain, and ``expand`` grows the neighboring domains to fill the gaps.

Filter options
    All :doc:`filter` options are accepted and applied to each tile and domain.

Branches
    ``--branch X`` (``-b``)
        Create a new branch: output results in ``{out}/{sample}/filterscan_{branch}``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`idmut` — produces the output this step consumes
- :doc:`clusterscan` — clusters the domains this step detects
- :doc:`filter` — run on each tile and domain internally
- :doc:`/formats/report/filterscan`
- :doc:`/use/inputs`, :doc:`/use/branch`, :doc:`/use/parallel`
