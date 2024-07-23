
Cluster: Infer alternative structures by clustering reads' mutations
--------------------------------------------------------------------------------

Cluster: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cluster input file: Mask report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of Mask report files as inputs for the Cluster step.
See :doc:`../inputs` for ways to list multiple files.

Cluster all masked reads in ``out``::

    seismic cluster out

Cluster reads from ``sample-1`` masked over reference reference ``ref-1``,
section ``abc``::

    seismic cluster out/sample-1/mask/ref-1/abc

Cluster: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _cluster_max:

Cluster setting: Maximum order (number of clusters)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To infer alternative RNA structures, SEISMIC-RNA uses an optimized version of
our original DREEM algorithm [`Tomezsko et al. (2020)`_], which is a type of
`expectation-maximization`_ (EM).
All EM algorithms need the order of clustering (i.e. number of clusters) to be
prespecified; however, the optimal order is unknown before the algorithm runs,
creating a `chicken-and-egg problem`_.

SEISMIC-RNA solves this problem by first running the EM algorithm at order 1,
then order 2, then 3, and so on until one of two limits is reached:

- The `Bayesian information criterion`_ (BIC) worsens upon increasing the order.
- The maximum order is reached.
  You can set this limit using ``--max-clusters`` (``-k``).
  If you run the entire workflow using ``seismic wf`` (see :doc:`./wf`), then
  the maximum order defaults to 0 (which disables clustering).
  If you run the Cluster step individually using ``seismic cluster``, then the
  maxmimum order defaults to 2 (the minimum non-trivial number).

.. note::
    If the BIC score worsens (increases) before reaching the maximum order,
    then clustering will stop.
    The report (see :doc:`../../formats/report/cluster`) records the maximum
    order you specified (field "Maximum Number of Clusters") and the order that
    yielded the best BIC (field "Optimal Number of Clusters"), which is always
    less than or equal to the maximum order you specified.

.. note::
    If you realize after clustering that it would have been better to have run
    clustering with a higher/lower maximum order, then you can edit the results
    using ``+addclust``/``+delclust`` (see :doc:`../adjclust`).

Cluster setting: Expectation-maximization iterations
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

`Expectation-maximization`_ is an iterative algorithm, meaning that it begins by
guessing an initial solution and then calculates progressively better solutions,
halting once successive solutions cease changing, which is called convergence.

You can limit the minimum/maximum number of iterations per number of clusters
using ``--min-em-iter`` and ``--max-em-iter``, respectively.
Generally, as the number of clusters increases, so does the number of iterations
required for convergence.
Thus, to treat different numbers of clusters more fairly, SEISMIC-RNA multiplies
the iteration limits by the number of clusters.
For example, if you use ``--max-em-iter 300``, then SEISMIC-RNA will allow up to
600 iterations for 2 clusters, 900 iterations for 3 clusters, and so on.
The exception is for 1 cluster: since all reads go into the same cluster, there
is no need to iterate, so the iteration limit is always the minimum possible, 2.

You can set the threshold for convergence with ``--em-thresh`` followed by the
minimum difference between log-likelihoods of successive iterations for the
iterations to be considered different.
For example, if you set the threshold to 0.1 with ``--em-thresh 0.1``, then if
iterations 38 and 39 had log-likelihoods of -7.28 and -7.17, respectively, then
the algorithm would keep going because their difference in log-likelihood (0.11)
would exceed the threshold; but if iteration 40 had a log-likelihood of -7.08,
then the algorithm would consider itself converged and stop running because the
difference in log-likelihood between iterations 40 and 39 would be 0.09, which
would be below the threshold.

Cluster setting: Expectation-maximization runs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

`Expectation-maximization`_ is guaranteed to return a locally optimal solution,
but there is no guarantee that the solution will be globally optimal.
To improve the odds of finding the global optimum, SEISMIC-RNA runs EM multiple
times (by default, 6 times), each time starting at a different initial guess.
The idea is that if multiple EM runs, initialized randomly, converge on the same
solution, then that solution is probably the global optimum.
You can set the number of independent EM runs using ``--em-runs`` (``-e``).

Cluster: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``OUT/SAMPLE/cluster/REFERENCE/SECTION``.

Cluster output file: Batch of cluster memberships
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Each batch of clustered reads contains a ``ClustBatchIO`` object and is saved to
the file ``cluster-batch-{num}.brickle``, where ``{num}`` is the batch number.
See :doc:`../../data/cluster/cluster` for details on the data structure.
See :doc:`../../formats/data/brickle` for more information on brickle files.

Cluster output file: Cluster report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA also writes a report file, ``cluster-report.json``, that records the
settings you used for running the Cluster step and summarizes the results, such
as the number of clusters, number of iterations, and the BIC scores.
See :doc:`../../formats/report/cluster` for more information.

.. note::
    You **must** look at the report file to determine whether your clusters come
    from true alternative structures or are just noise and artifacts.
    See :ref:`clust_verify` for how to verify that your clusters are real.

.. _clust_verify:

Cluster: Verify clusters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You **must** check whether your clusters are real or artifacts.

In your cluster report:

- The number of clusters that SEISMIC-RNA found is Optimal Number of Clusters.
  Several important caveats exist about this number:

  - This number can never exceed the Maximum Number of Clusters.
    So if you want to know whether an RNA forms *N* alternative structures, the
    results of clustering can provide useful information only if you set the
    Maximum Number of Clusters to at least *N*.
  - A "cluster" is as subjective as a "conformational state": two clusters can
    correspond to completely different structures at one extreme and to slightly
    different structures at the other.
    With more reads comes better ability to distinguish clusters that are more
    similar -- the same way that, in a study examining differences between two
    groups, larger sample sizes would enable finding more subtle differences.
    Thus, the number of clusters you find will generally increase with more
    reads, but that doesn't mean that your RNA actually forms more structures,
    just that you can resolve more subtle structural differences.
  - The Number of Unique Bit Vectors is the number of reads that were used for
    clustering; it should be about 20,000 at minimum, and ideally ≥ 30,000.
    If you have < 20,000 unique bit vectors, then clustering will probably not
    be able to find real clusters; so if the Optimal Number of Clusters is 1,
    then that does not mean your RNA necessarily forms only one structure.

- `Expectation-maximization`_ is guaranteed to find a local optimum, but not a
  global optimum.
  SEISMIC-RNA thus runs multiple trajectories from different starting points; if
  the trajectories converge to the same solution, then that solution is likely
  (but still not necessarily) the global optimum.
  You must check if your trajectories converged to the same solution by checking
  the fields "NRMSD from Run 0" and "Correlation with Run 0" in the report.
  If all runs converged to identical solutions, then every NRMSD would be 0 and
  every Correlation would be 1.
  Generally, the runs are sufficiently reproducible if runs 1 and 2 have NRMSDs
  less than 0.05 and Correlations greater than 0.98 with respect to run 0.
  If not, then there you have no evidence that run 0 is the global optimum for
  that number of clusters, so it would be best to rerun clustering using more
  independent runs to increase the chances of finding the global optimum.

Cluster: Troubleshoot and optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run Cluster with higher orders, without repeating the work already done
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The tool ``+addclust`` exists for this purpose: see :ref:`addclust`.

Delete unnecessary higher orders, without repeating the work already done
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The tool ``+delclust`` exists for this purpose: see :ref:`delclust`.

Cluster takes too long to finish
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Adjust the settings of ``seismic cluster``:

  - Increase the threshold for convergence (``--em-thresh``).
    Larger thresholds will make clustering converge in fewer iterations at the
    cost of making the runs end at more variable solutions.
    Check the Log Likelihood per Run field to verify that clustering is finding
    the global optimum; see :ref:`clust_verify` for more information.
  - Decrease the number of independent runs (``--em-runs``/``-e``) to 3 or 4;
    don't go below 3 for anything you intend to publish, or else you won't be
    able to tell if your clustering is finding the global optimum.

.. _Tomezsko et al. (2020): https://doi.org/10.1038/s41586-020-2253-5
.. _expectation-maximization: https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm
.. _chicken-and-egg problem: https://en.wikipedia.org/wiki/Chicken_or_the_egg
.. _Bayesian information criterion: https://en.wikipedia.org/wiki/Bayesian_information_criterion
