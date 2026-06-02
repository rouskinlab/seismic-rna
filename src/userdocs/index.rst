########################################################################
SEISMIC-RNA
########################################################################

Structure Ensemble Inference by Sequencing, Mutation Identification, and
Clustering of RNA


News
================================================================================

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
