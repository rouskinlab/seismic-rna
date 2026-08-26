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
Optionally (``--validate-gaps``), it then checks whether each pair of
neighboring domains is truly independent or was mistakenly split apart by
filterscan, merging and re-clustering any domains that turn out to be coupled.
Use it after :doc:`filterscan` to resolve the structural heterogeneity within
each data-driven domain.


How it works
================================================================================

Checking whether adjacent domains are really separate (optional)
    filterscan cuts a domain boundary wherever bridge pairs (the
    anti-correlated pairs of positions it uses as evidence of alternative
    structures) become locally sparse. Occasionally, that depletion happens
    in the middle of what is really a single domain -- for instance, because
    that stretch happens to have fewer testable pairs -- so filterscan
    reports one true domain as two separate, adjacent domains. With
    ``--validate-gaps`` turned on, for every pair of adjacent domains that
    each resolved into 2 or more clusters, clusterscan looks at the reads
    that span the gap between them (reads that reach into both domains at
    once) and tests whether the alternative structures found on each side
    occur together about as often as expected if the two sides fold
    independently, or whether certain combinations turn up unexpectedly
    often or rarely together -- a sign that the "two" domains are really one
    coupled unit that got cut in half. ``--gap-min-assoc`` sets how much
    apparent coupling is needed before a gap is judged spurious: higher
    values merge only strongly coupled domains; lower values merge more
    readily. Gaps with too few spanning reads to test reliably, or where
    either side settled on just 1 cluster (so there is nothing to compare),
    are left alone.

    When a gap looks spurious, clusterscan merges the two domains into one
    and re-clusters the merged span from scratch, since the domain's true
    clusters likely correspond to particular pairings of structure across
    the whole merged region, not to the halves it was mistakenly cut into.
    Because merging can change what looks coupled at neighboring gaps,
    clusterscan then re-checks every remaining gap and, among any that still
    look spurious, merges the most strongly coupled pair first; it repeats
    this until every gap between the remaining domains looks like a genuine,
    independent break. ``--validate-gaps`` is off by default because it can
    be expensive: each merge re-filters and re-clusters the combined region.
    The final report records which of filterscan's original domains were
    merged into each final domain; see :doc:`/formats/report/clusterscan`.


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
    Summary of the run, including the best number of clusters for each final
    domain and, if ``--validate-gaps`` merged any domains, which of
    filterscan's original domains were combined into each one.
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

Domain merging
    ``--validate-gaps/--no-validate-gaps``
        Check every gap between adjacent domains for whether the two domains
        are actually coupled, merging and re-clustering them if so (default:
        off, since it can be expensive). See "How it works" above.
    ``--gap-min-assoc F``
        Minimum apparent coupling between two domains' structures required
        to judge a gap spurious and merge them (default 0.1). Higher values
        merge only strongly coupled domains; lower values merge more readily.

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
