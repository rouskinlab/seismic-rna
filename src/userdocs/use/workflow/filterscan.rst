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


How it works
================================================================================

Finding correlated pairs
    For every pair of positions close enough together (within
    ``--band-width``, if set) and covered together by enough reads
    (``--min-pair-coverage``), SEISMIC-RNA checks whether the two positions
    mutate together more or less often than would be expected if they mutated
    independently. Positions in regions that form multiple structures tend to
    show positive correlations if they are both paired in one structure and both
    unpaired in another, and negative correlations if each is unpaired in a
    different structure. Two postitions located in a region that forms only one
    structure or located in different independently-folding structures are
    expected to show no correlation. Each pair of positions is scored by how
    strongly its co-mutation pattern exceeds what independence would predict.

Grouping pairs into domains
    A domain is a stretch of the RNA where many pairs of positions show
    elevated correlation with each other. SEISMIC-RNA searches for the way
    of dividing the RNA into domains (and background, i.e. not part of any
    domain) that best explains the observed pattern of correlated pairs, so
    that positions within the same domain tend to be correlated with each
    other while positions in different domains do not.

Telling real domains from chance
    Any dataset shows some correlation between positions purely by chance,
    even with no real structure, so SEISMIC-RNA needs a way to judge how
    much correlation is too much to be coincidence. For every candidate
    domain, it computes an exact statistical test of whether the observed
    correlation could plausibly have arisen if every position in the
    candidate mutated independently, and corrects for the fact that it is
    testing a great many overlapping candidate domains at once (so that a
    domain is called only when it is a genuinely unusual result, not merely
    the best of many chance draws). ``--detect-fdr`` sets the target false
    discovery rate (FDR) for this initial call: how willing SEISMIC-RNA is
    to call a domain. It intentionally has an unusually high default for an
    FDR (0.1) to make it more sensitive because here, false negatives are
    worse than false positives. A false positive (a detected domain that
    really forms only one structure) merely slows down the workflow:
    ClusterScan needs to spend time clustering it, but a false positive
    domain would likely yield 1 cluster due to the Cluster step's stringent
    filters, so the final result will likely be correct. A false negative
    (failing to detect where the RNA really forms multiple structures)
    would not be passed into the Cluster step at all and hence the final
    result would incorrectly be 1 cluster.

Merging domains across gaps
    Sometimes two domains are separated by a short, unstructured stretch
    (for example, an unpaired linker) that shows no correlation of its own,
    yet the two domains are nonetheless connected by real, direct
    long-range pairs (for example, a helix whose two strands lie on either
    side of the gap). SEISMIC-RNA checks, at every point within the gap,
    whether pairs of positions spanning that specific point still show a
    real, direct correlation, judged at its own target false discovery rate,
    ``--merge-fdr``; if the connection holds all the way across the gap, it
    merges the two domains into one. This keeps a single, genuinely long-range
    structure from being cut into pieces just because an unstructured stretch
    happens to separate the two ends where its correlated pairs anchor, while
    still splitting apart two domains that are truly independent, with nothing
    but coincidental noise between them.


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
    ``--min-expect-both N``
        Analyze only pairs of positions whose expected number of jointly
        mutated reads, if the two positions mutated independently, is at
        least this value (default 5): standard practice for the
        statistical test SEISMIC-RNA uses, which becomes unreliable when
        this expected count drops too low.
    ``--anticorr-only/--no-anticorr-only``
        Count only pairs of positions that are negatively correlated, i.e.
        mutated less often together than expected by chance (default: only
        negatively correlated pairs). Two positions are negatively
        correlated when they belong to different alternative structures: a
        read that has one modified rarely has the other modified too, since
        each read reflects only one structure at a time. Positively
        correlated pairs (mutated together more often than chance) can
        arise from other sources, such as local structural or chemical
        effects unrelated to alternative structures, so they are excluded
        by default. Disable this option to consider both kinds of
        correlation, as SEISMIC-RNA did in earlier versions.
    ``--detect-fdr F``
        How willing to be to call a region a domain, expressed as a false
        discovery rate (default 0.1): SEISMIC-RNA corrects for the fact
        that it tests a great many overlapping candidate domains for
        correlation exceeding chance. Higher values call more (and weaker)
        domains; lower values call fewer, more conservative domains.
    ``--merge-fdr F``
        How willing to be to join two domains separated by a gap, expressed
        as a false discovery rate (default 0.1): SEISMIC-RNA corrects for
        the fact that it tests every point within the gap for a direct
        correlation crossing it. Higher values merge more (and weaker)
        connections; lower values merge fewer, more conservative ones.

Domain length filters
    ``--min-domain-length N``
        Keep only the domains with at least this many positions (default 20).

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
