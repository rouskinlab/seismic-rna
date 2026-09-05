********************************************************************************
seismic filterscan
********************************************************************************


Purpose
================================================================================

``seismic filterscan`` searches a full-length RNA for **domains**: regions that
appear to fold independently, revealed by pairs of positions whose mutations
are anti-correlated with each other (bridge pairs).
It slides overlapping tiles along the RNA, runs the filter step on each tile,
finds pairs of positions that mutate together less often than expected by
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
    are mutated together *less* often than would be expected if they
    mutated independently, i.e. whether they are anti-correlated. Two
    positions that belong to different alternative structures tend to be
    anti-correlated: a read that has one of them modified rarely has the
    other modified too, since any one read reflects only one structure at a
    time. Only anti-correlation is used as evidence of a domain boundary;
    the opposite pattern, positive correlation, is not, because it can also
    arise from causes that have nothing to do with alternative structures,
    such as reads with more RNA overall picking up more mutations at every
    position. A pair that clears this bar becomes a **bridge pair** only if
    its anti-correlation is also statistically significant, correcting for
    the very large number of pairs being tested at once (``--pair-fdr``),
    and large enough in size to matter biologically (``--min-fold-change``):
    a real split between alternative structures changes how often both
    positions are modified in the same read by a specific amount that
    depends on the split (for instance, a 50/50 split between two
    structures depletes joint modifications differently than a 90/10
    split), so this second check screens out pairs that are "significant"
    only because they were measured with enormous read depth, not because
    the effect is actually large.

Grouping pairs into domains
    A domain is a stretch of the RNA with an unusually high concentration
    of bridge pairs. SEISMIC-RNA searches for the single way of dividing
    the RNA into domains (and background, i.e. not part of any domain) that
    best explains the observed bridge pairs overall, rather than greedily
    accepting the first plausible-looking domain it finds. Each candidate
    domain is judged only against its own bridge pairs, so a strongly
    structured domain in one part of the RNA cannot make an unrelated,
    distant stretch look domain-like too (or hide a real domain next to it).

Telling real domains from chance
    Even outside of any real domain, a few pairs will look like bridge
    pairs purely by chance, so SEISMIC-RNA needs a way to judge how many
    bridge pairs in a candidate domain is too many to be coincidence. It
    first estimates the background rate at which bridge pairs turn up
    outside of any real domain, then, for every candidate domain, computes
    an exact statistical test of whether its own concentration of bridge
    pairs is higher than that background rate would predict, correcting
    for the fact that it is testing a great many overlapping candidate
    domains at once (so that a domain is called only when it is a
    genuinely unusual result, not merely the best of many chance draws).
    ``--detect-fdr`` sets the target false discovery rate (FDR) for this
    initial call: how willing SEISMIC-RNA is to call a domain. It
    intentionally has an unusually high default for an FDR (0.1) to make it
    more sensitive because here, false negatives are worse than false
    positives. A false positive (a detected domain that really forms only
    one structure) merely slows down the workflow: ClusterScan needs to
    spend time clustering it, but a false positive domain would likely
    yield 1 cluster due to the Cluster step's stringent filters, so the
    final result will likely be correct. A false negative (failing to
    detect where the RNA really forms multiple structures) would not be
    passed into the Cluster step at all and hence the final result would
    incorrectly be 1 cluster.

Merging domains across gaps
    Sometimes two domains are separated by a short, unstructured stretch
    (for example, an unpaired linker) that has too few bridge pairs of its
    own to look like part of either domain, yet the two domains are
    nonetheless connected by real, direct long-range bridge pairs (for
    example, a helix whose two strands lie on either side of the gap).
    SEISMIC-RNA checks, at every point within the gap, whether the bridge
    pairs spanning that specific point are still enough to look like a real
    connection rather than chance noise, judged at its own target false
    discovery rate, ``--merge-fdr``; if the connection holds all the way
    across the gap, it merges the two domains into one. This keeps a
    single, genuinely long-range structure from being cut into pieces just
    because an unstructured stretch happens to separate the two ends where
    its bridge pairs anchor, while still splitting apart two domains that
    are truly independent, with nothing but coincidental noise between them.

Growing domain edges to catch connections just outside them
    Even after the steps above settle on a domain's boundaries, a domain
    can still have a handful of real bridge pairs reaching in from just
    past one of its edges -- for example, a long-range helix whose far
    strand lies just outside the boundary, connected to the rest of the
    domain across a short stretch that, on its own, is too sparse in bridge
    pairs to have been swept into the domain by the earlier steps. Checked
    one at a time, whether any single one of these edge-anchored pairs
    looks like a real connection is unstable: it can flip from "yes" to
    "no" from one position to the next simply by chance. Checked together
    as a group, however, a genuine cluster of them is unambiguous. So
    SEISMIC-RNA looks at the farthest bridge pair reaching in from outside
    a domain's edge and asks, using the same statistical test as
    ``--merge-fdr`` above, whether everything reaching in from that
    position to the domain, taken as a whole, is a real connection rather
    than chance noise. If so, it grows the domain out to that pair and
    repeats the check, since growing one edge can expose further bridge
    pairs just beyond it; if the farthest pair turns out to be an isolated
    coincidence with nothing else connecting it to the domain, the edge is
    left where it was. A domain grown this way never extends past
    ``--max-domain-length`` or into a neighboring domain.


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
    The bridge pairs of positions found in the RNA, with each pair's P-value,
    BH-adjusted P-value, and fold change.

``domains.csv``
    The coordinates (5' and 3' ends) of the final domains, and how each was
    produced: called as-is, grown by ``--widen``, or inserted by ``--fill``.

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
    ``--pair-fdr F``
        How willing to be to count an anti-correlated pair as a bridge
        pair, expressed as a false discovery rate (default 0.05):
        SEISMIC-RNA corrects for the fact that it tests every eligible pair
        of positions at once. Higher values call more (and weaker) bridge
        pairs; lower values call fewer, more conservative ones.
    ``--min-fold-change F``
        Minimum fold change, between the number of reads expected to be
        modified at both positions of a pair (if the two positions were
        independent) and the number actually observed, for the pair to
        count as a bridge pair (default 2.0). Screens out pairs whose
        anti-correlation is statistically significant only because they
        were measured with very high read depth, not because the effect on
        the underlying structures is large. Higher values require a larger
        effect.
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
    ``--max-domain-length N``
        Cap the length of every domain, in positions (default 0 = 2× the
        median read length): the search for domains never calls one longer
        than this, and ``--widen``/``--fill`` (below) never grow or insert
        one longer than this either, so every domain fits within whatever
        window the Cluster step will later analyze.
    ``--min-domain-length N``
        Keep only the domains with at least this many positions (default
        20): drops the shortest, least confident calls. Skipped when
        ``--widen`` is set, since widening grows a short domain into its
        neighboring gap instead of dropping it.

Gap handling
    ``--widen/--no-widen``
        Grow each domain into the unstructured gaps on either side of it
        (default: off), up to ``--max-domain-length`` in total. If two
        neighboring domains both reach that limit before their shared gap
        closes, whatever space is left between them stays a gap (which
        ``--fill``, below, can then close).
    ``--fill/--no-fill``
        Insert domains into every gap that's left -- including the very
        ends of the scanned region -- so every position ends up belonging
        to exactly one domain (default: off). A gap longer than
        ``--max-domain-length`` is split into the smallest number of
        equally-sized domains that each stay within that limit, rather than
        becoming one long, unwieldy domain. If ``--widen`` is also set, it
        runs first, and ``--fill`` closes whatever gaps it leaves behind.

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
