########################################################################
SEISMIC-RNA
########################################################################

Structure Ensemble Inference by Sequencing, Mutation Identification, and
Clustering of RNA


News
================================================================================

.. admonition:: Version 0.26 released
   :class: news

   SEISMIC-RNA 0.26 adds cofolding of two RNAs, replaces the ``ensembles``
   step with two separate steps, and shows progress bars while it runs.
   Output directories from version 0.25 can be read by version 0.26 without
   migrating, except for output of the ``ensembles`` step, which must be
   regenerated with the two steps that replace it.

   **New features**

   - :doc:`seismic duplex </use/workflow/duplex>` combines two references
     into one duplex, and :doc:`seismic fold </use/workflow/fold>` then
     predicts the structure the two molecules form together (cofolding with
     RNAcofold).  A second strand can be another dataset, the same dataset
     (``--dimer``, for a homodimer), or a bare sequence with no data
     (``--duplex-file`` / ``--duplex-sequence``), such as an antisense
     oligonucleotide.
   - Progress bars now show what SEISMIC-RNA is working on, one line per
     task, for every parallelized step and for :doc:`seismic wf
     </use/workflow/wf>` as a whole.  Hide them with ``--no-progress``
     (they are hidden automatically with ``-q`` or when stderr is not a
     terminal).
   - :doc:`seismic collate </use/workflow/collate>` has been rewritten and
     now runs by default at the end of :doc:`seismic wf </use/workflow/wf>`.
   - :doc:`seismic filterscan </use/workflow/filterscan>` detects domains
     with a new model based on bridge pairs — the position pairs that link
     two parts of a domain — tuned with ``--pair-fdr``,
     ``--min-fold-change``, ``--detect-fdr``, ``--merge-fdr``, and
     ``--min-pairs``.
   - SEISMIC-RNA starts up faster, and its log messages have been simplified
     from eight levels to five.

   **Changes to be aware of**

   - The ``ensembles`` step has been split into two steps:
     :doc:`filterscan </use/workflow/filterscan>` (detect domains) and
     :doc:`clusterscan </use/workflow/clusterscan>` (cluster each domain).
     Its output directory ``ensembles/`` is likewise replaced by
     ``filterscan/`` and ``clusterscan/``.
   - Graphs are now written only as interactive HTML files.
     The ``--svg``, ``--pdf``, and ``--png`` options have been removed, along
     with the Kaleido dependency they required; ``--html`` / ``--no-html``
     controls whether each graph is written.
   - The domain-finding options of the old ``ensembles`` step
     (``--threshold-divisor``, ``--min-cluster-length``,
     ``--max-cluster-length``, and ``--gap-mode``) have been removed in
     favor of the new options of
     :doc:`filterscan </use/workflow/filterscan>`, and ``--min-pairs`` now
     sets the minimum number of bridge pairs per domain.
   - SEISMIC-RNA now requires Python 3.13 (previously 3.11), and its
     dependencies have been updated to match.

   See the `changelog`_ for the complete list of changes.

.. admonition:: Version 0.25 released — breaking changes require migration
   :class: news

   SEISMIC-RNA 0.25 contains several backwards-incompatible changes.
   Existing output directories produced by version 0.24 cannot be read by
   version 0.25 without first running the migration command described below.

   **What changed**

   - The ``relate`` step has been renamed to **IDmut** (``idmut``).
     Output subdirectories previously named ``relate/`` are now ``idmut/``;
     file names that began with ``relate-`` now begin with ``idmut-``.
   - The ``mask`` step has been renamed to **Filter** (``filter``).
     Output subdirectories previously named ``mask/`` are now ``filter/``;
     file names that began with ``mask-`` now begin with ``filter-``.
   - Several command-line options have been renamed or removed.
     Report JSON fields have been updated to match.

   **How to migrate**

   Use :doc:`seismic migrate </use/utility/migrate>` to update an existing
   output directory to the version 0.25 format::

       seismic migrate out -o out-new

   where ``out`` is your old output directory (from version 0.24) and
   ``out-new`` is a new directory that will be created to hold the updated
   outputs.  The original ``out`` directory is never modified; if an error
   occurs the incomplete ``out-new`` directory is removed automatically.

   Once the command succeeds, verify your results using ``out-new`` in
   place of ``out`` for any downstream steps.


.. toctree::
    :maxdepth: 1

    why/index
    works/index
    install/index
    use/index
    tutorials/index
    issues
    about/index
    cli
    api/index
    formats/index
    data/index
    algos/index
    faqs/index

.. _changelog: https://github.com/rouskinlab/seismic-rna/blob/main/CHANGELOG.md
