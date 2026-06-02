********************************************************************************
seismic fold
********************************************************************************


Purpose
================================================================================

``seismic fold`` uses per-position mutation rates to predict RNA secondary
structures using one of two external folding programs.
It normalizes the mutation rates into reactivity scores and passes them as
soft constraints to guide the structure prediction.

Requires either RNAstructure_ (Fold/ShapeKnots) or ViennaRNA_ (RNAfold/RNAsubopt)
depending on the ``--fold-backend`` setting.


Inputs
================================================================================

Table CSV files or table output directories
    Per-position mutation rate tables produced by ``seismic table``.
    Pass a table CSV, the directory containing it, or a higher-level directory.
    See :doc:`/use/inputs`.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/fold/{ref}/{reg}/``.

``{profile}.ct``
    One CT file per mutational profile, named ``{reg}__average.ct`` (or
    ``{reg}__cluster-{K}-{k}.ct`` for clustered data); contains the predicted
    base pairs. See :doc:`/formats/data/ct`.

``fold-report.json``
    Summary of settings and results.
    See :doc:`/formats/report/fold`.


Quick example
================================================================================

Fold filter output using the default DMS settings::

    seismic fold out/sample-1/filter/ref-1


Options
================================================================================

Backend
    ``--fold-backend {auto|rnastructure|viennarna}``
        Folding program to use.
        ``auto`` (default): uses RNAstructure for DMS, ViennaRNA otherwise.
    ``--pseudoknots/--no-pseudoknots``
        Predict pseudoknotted structures using ShapeKnots (requires
        RNAstructure; default off).

Reactivity normalization
    ``--fold-quantile F``
        Winsorize reactivities to this quantile before folding (default 0.95).
    ``--fold-energy-method {auto|cordero|deigan|eddy}``
        How to convert reactivities to folding energy penalties (default auto).
    ``--deigan-slope F``
        Slope (kcal/mol) for the Deigan energy method (default 1.8).
    ``--deigan-intercept F``
        Intercept (kcal/mol) for the Deigan energy method (default −0.6).

Folding conditions
    ``--fold-temp F``
        Temperature in Celsius (default 37.0).
    ``--fold-md N``
        Maximum base-pair distance in nucleotides; 0 for no limit (default 0).
    ``--fold-isolated/--fold-stacked``
        Allow isolated (non-stacked) base pairs (default off).

Output structures
    ``--fold-mfe/--fold-sub``
        Predict only the minimum-free-energy structure (default off = predict
        suboptimal structures as well).
    ``--fold-max N``
        Maximum number of suboptimal structures to output (default 20).
    ``--fold-percent F``
        Stop outputting structures when energy exceeds MFE by this % (default 20).

Constraints
    ``--fold-constraint FILE``
        Force specific bases to be paired or unpaired from a constraint file.

Region selection
    ``--fold-coords REF FIRST LAST``
        Fold only positions FIRST–LAST of REF.
    ``--fold-primers REF FWD REV``
        Define the region by primer sequences.
    ``--fold-regions-file FILE``
        Define regions to fold from a CSV file.

Branches
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/fold_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

All structures have very few base pairs
    Mutation rates may be too high or not normalized correctly.
    Check that the probe preset (DMS vs. SHAPE) and ``--fold-quantile`` are
    appropriate for your data.

Fold fails immediately
    The required backend (RNAstructure or ViennaRNA) may not be installed.
    Check that ``Fold``, ``ShapeKnots``, ``RNAfold``, or ``RNAsubopt`` is on
    your PATH.


See also
================================================================================

- :doc:`table` — produces the mutation rate tables this step uses
- :doc:`draw` — visualizes the predicted structures
- :doc:`/formats/data/ct`, :doc:`/formats/report/fold`
- :doc:`/use/branch`, :doc:`/use/parallel`


.. _RNAstructure: https://rna.urmc.rochester.edu/RNAstructure.html
.. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/
