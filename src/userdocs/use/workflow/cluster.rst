********************************************************************************
seismic cluster
********************************************************************************


Purpose
================================================================================

``seismic cluster`` fits a mixture model to the filtered reads, finding
subpopulations (clusters) with distinct mutation patterns.
Use it when you expect your RNA preparation to contain multiple structural
conformations or when DMS-MaPseq reveals heterogeneous reactivity.


Inputs
================================================================================

Filter output directories or report files
    One or more paths to ``seismic filter`` output.
    See :doc:`/use/inputs`.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/cluster/{ref}/{reg}/``.

``cluster-batch-{num}.brickle``
    Per-cluster read assignments for one batch.
    See :doc:`/formats/data/brickle`.

``cluster-report.json``
    Summary of settings and the best number of clusters found.
    See :doc:`/formats/report/cluster`.


Quick example
================================================================================

Cluster filter output, searching up to 3 clusters::

    seismic cluster -k 3 out/sample-1/filter/ref-1


Options
================================================================================

Number of clusters
    ``--max-clusters N`` (``-k``)
        Maximum number of clusters to try (default 0 = no limit; set this).
    ``--min-clusters N``
        Minimum number of clusters to try (default 1).
    ``--try-all-ks/--stop-best-k``
        Try all values of K even after finding the apparent best (default off).
    ``--write-all-ks/--write-best-k``
        Write results for all K values, not just the best (default off).

EM algorithm
    ``--min-em-runs N`` (``-e``)
        Minimum successful EM runs per K (default 6).
    ``--max-em-runs N`` (``-E``)
        Maximum EM attempts per K (default 30).
    ``--min-em-iter N``
        Minimum EM iterations per run (default 10).
    ``--max-em-iter N``
        Maximum EM iterations per run (default 500).
    ``--em-thresh F``
        Stop EM when log-likelihood improvement falls below F (default 0.37).

Run quality filters (used to discard low-quality EM solutions)
    ``--max-pearson-run F``
        Discard runs where two clusters correlate above F (default 0.9).
    ``--min-marcd-run F``
        Discard runs where two clusters are less than F apart (default 0.016).
    ``--max-gini-run F``
        Discard runs where any cluster has a Gini coefficient above F
        (default 0.667).

K selection filters
    ``--min-pearson-vs-best F``
        Remove K values where every run correlates below F with the best K
        (default 0.97).
    ``--max-marcd-vs-best F``
        Remove K values where every run is more than F from the best K
        (default 0.008).

Jackpotting
    ``--jackpot/--no-jackpot``
        Check for over-represented reads caused by PCR jackpotting (default on).
    ``--max-jackpot-quotient F``
        Discard runs whose jackpotting quotient exceeds F (default 1.1).

Branches
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/cluster_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

Clustering always finds 1 cluster
    The data may not have enough reads or positional coverage to separate
    clusters.
    Check that ``--min-ninfo-pos`` in filter is not too strict, or that you
    have enough reads.

Clustering is very slow
    Reduce ``--max-em-runs`` or ``--max-em-iter``.
    Also reduce ``--num-cpus`` if memory is the bottleneck.


See also
================================================================================

- :doc:`filter` — produces the data this step consumes
- :doc:`join` — combine cluster results across regions
- :doc:`/formats/report/cluster`
- :doc:`/use/parallel`, :doc:`/use/branch`
