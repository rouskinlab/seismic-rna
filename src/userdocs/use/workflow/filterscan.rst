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
    ``--band-width N``
        Consider only pairs of positions no farther apart than this many
        bases when looking for domains (default 0 = no extra limit beyond
        the tile length).
    ``--min-pair-coverage N``
        Analyze only pairs of positions with at least this many jointly
        covering reads (default 1000): pairs with less coverage are too
        noisy to score reliably.
    ``--domain-fdr F``
        How willing to be to call a region a domain, expressed as a false
        discovery rate (default 0.1): SEISMIC-RNA simulates data with no
        real structure and compares it with the real data to judge how
        much correlation could arise by chance alone. Higher values call
        more (and weaker) domains; lower values call fewer, more
        conservative domains.
    ``--n-null-replicates N``
        Number of simulated no-structure datasets to use for calibrating
        ``--domain-fdr`` (default 10).
    ``--seed N``
        Seed for the random number generator used to simulate the
        no-structure datasets (default: none, i.e. nondeterministic). Set
        a value to make domain calling reproducible.

Domain length filters
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
