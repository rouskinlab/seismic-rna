
Fold: Predict RNA secondary structures using mutation rates
--------------------------------------------------------------------------------

Fold: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fold input file: Filter or Cluster positional table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of positional table files of masked or clustered reads
(``filter-per-pos.csv`` or ``clust-per-pos.csv``, respectively) as inputs.
See :doc:`../inputs` for ways to list multiple files.
(SEISMIC-RNA will not crash if you give other type of table files, such as a
``idmut-per-pos.csv`` or ``filter-per-read.csv.gz`` file, but will ignore them.)

To predict structures using the mutational profiles in all valid tables in the
directory ``{out}``, you could use the command ::

    seismic fold {out}

Fold: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fold setting: Choose a folding backend
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``seismic fold`` supports three folding backends, selected with
``--fold-backend``:

=================== =================== ========================================
``--fold-backend``  Program             Package
=================== =================== ========================================
``Fold`` (default)  RNAstructure Fold   RNAstructure_ (≥ 6.6)
``ShapeKnots``      RNAstructure        RNAstructure_ (≥ 6.6) — predicts
                    ShapeKnots          pseudoknots
``RNAFold``         ViennaRNA RNAfold   ViennaRNA_ (≥ 2.7.2) — see
                                        :ref:`install_dependencies`
=================== =================== ========================================

All three backends accept normalized mutation rates as soft constraints to guide
structure prediction (see `Fold setting: Energy method`_).
ShapeKnots is the only backend that can predict pseudoknots.
RNAFold is the only backend that supports the Eddy energy method.

Fold setting: Energy method
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Mutation rates are incorporated into the folding energy function as soft
constraints.
Use ``--fold-energy-method`` to choose the method:

- ``Deigan`` (default): SHAPE pseudo-energy term ``m * log(reactivity + 1) + b``
  with slope ``--deigan-slope`` (default 1.8 kcal/mol) and intercept
  ``--deigan-intercept`` (default −0.6 kcal/mol).
  Works with all three backends: RNAstructure uses SHAPE-directed folding;
  ViennaRNA uses a soft-constraint file.

- ``Cordero``: Hard partition of positions into paired/unpaired constraints
  based on a reactivity threshold.
  Requires ``--fold-backend Fold`` or ``--fold-backend ShapeKnots``.

- ``Eddy``: Uses ViennaRNA's built-in soft constraint facility.
  Requires ``--fold-backend RNAFold``.

Fold setting: Define regions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can predict structures of the full reference sequences or specific regions.
See :doc:`../regions` for ways to define regions.

Defining regions in ``seismic fold`` works identically to ``seismic filter`` but
accomplishes a very different purpose.
Regions in ``seismic fold`` determine for which parts of the reference sequence
to predict structures.
Regions in ``seismic filter`` determine for which parts of the reference sequence
to use mutational data.
SEISMIC-RNA allows these regions to be different.
There are several common scenarios:

- The region you are folding matches the region for which you have data.
  For example, you could have mutationally profiled a full transcript and now
  want to predict the structure of the full transcript using the data from the
  full mutational profile.
- You are folding a region that contains and is longer than the region for
  which you have data.
  For example, you could have mutationally profiled a short amplicon from a much
  longer transcript; and after clustering that amplicon, you want to model each
  alternative structure of the long transcript while using the short mutational
  profile of each cluster to guide the structure predictions.
- You are folding a short region that is contained by a longer region for
  which you have mutational profiling data.
  For example, you could have mutationally profiled a full transcript and now
  want to predict the structure of a small part of the transcript that you are
  reasonably sure does not interact with any other part of the transcript.

Fold setting: Quantile for normalization
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Folding requires that the mutation rates be normalized to the interval [0, 1].
Use ``--fold-quantile`` (default 0.95) to set the quantile to which reactivities
are normalized and winsorized before folding.
See :doc:`../normalize` for more information on normalization.

Fold setting: RNAstructure parameters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``seismic fold`` exposes several options for the RNAstructure Fold and
ShapeKnots programs (see the `documentation for Fold`_ for details on each
option).
Options marked with (†) are also honoured by the ViennaRNA RNAfold backend.

============================== =========================== ===========================================================================================
Option in ``seismic fold``     Option in RNAstructure      Brief explanation
============================== =========================== ===========================================================================================
``--fold-temp``                ``--temperature``           temperature (°C) of folding (default 37)
``--fold-constraint``          ``--constraint``            optional `folding constraints file`_ with forced pairs/unpaired bases (†)
``--fold-md`` (†)              ``--maxdistance``           maximum distance between paired bases (0 = no limit)
``--fold-mfe`` (†)             ``--MFE``                   predict only the optimal structure (same result as ``--fold-max 1``, but about twice as fast)
``--fold-max`` (†)             ``--maximum``               maximum number of structures to predict (ignored if using ``--fold-mfe``)
``--fold-percent``             ``--percent``               maximum % difference in free energy of predicted structures (ignored if using ``--fold-mfe``)
``--fold-isolated``            ``--maxloop``-adjacent      allow isolated (non-stacked) base pairs (default: disallowed)
============================== =========================== ===========================================================================================

Fold setting: ViennaRNA-specific parameters
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When ``--fold-backend RNAFold`` is selected, you may additionally pass a
commands file to RNAfold using ``--fold-commands``.
This file is forwarded verbatim to RNAfold's ``--commands`` option.

Fold setting: Dry run
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Use ``--fold-dry-run`` to generate the input files and command for the folding
backend without actually running it.
This is useful for inspecting the soft-constraint data files or debugging the
command line before a long production run.
Use ``--fold-real-run`` (the default) to run folding normally.

Fold: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{sample}/fold/{ref}/{reg}``,
where ``{out}`` is the output directory, ``{sample}`` is the sample, ``{ref}``
is the reference, and ``{reg}`` is the region you *folded* (**not** that from
which the data came).
The files for each predicted structure are named ``{reg}__{profile}``, where
``{reg}`` is the region from which the *data* came (**not** that which you
folded) and ``{profile}`` is the mutational profile of those data, which can be
``average`` (ensemble average) or ``cluster-{n}-{i}`` (where ``{n}`` is the
number of clusters and ``{i}`` is the cluster number).

Fold output file: Fold report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA writes a report file, ``fold-report.json``, to record the settings
you used for running the Fold step.
See :doc:`../../formats/report/fold` for more information.

Fold output file: Connectivity table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The primary output is a connectivity table file.
For details on this format, see :doc:`../../formats/data/ct`.

.. _fold_db:

Fold output file: Dot-bracket structure
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Fold step also outputs the structures in dot-bracket format, which you can
copy-paste into RNA drawing software such as `VARNA`_.
For details on this format, see :doc:`../../formats/data/db`.

Fold output file: VARNA color file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Fold step also outputs the normalized mutation rates in VARNA color format,
which you can import into the RNA drawing software `VARNA`_.
For details on this format, see :doc:`../../formats/data/varna-color`.

Fold: Visualize structures in VARNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`VARNA`_ is a third-party application for drawing RNA structures.
To draw a structure from SEISMIC-RNA in VARNA:

1.  Install (if needed) and launch VARNA.
2.  Open your dot-bracket file (see :ref:`fold_db`) in a text editor.
3.  Right-click the drawing canvas, select "File" > "New...", and copy-paste the
    sequence and dot-bracket structure.
4.  Adjust the layout of the structure by clicking and dragging.
5.  To color the bases by their mutation rates, right-click the drawing canvas,
    select "Display" > "Color map" > "Load values...", copy-paste the path to
    your VARNA color file into the box or click "Choose file" and navigate to
    your VARNA color file, and click "OK" to load the file.
6.  To customize the colors, select "Display" > "Color map" > "Style...":

    - Drag a color bar to adjust its location.
    - Click the square below a color bar to change its color.
    - Click the X below the square to delete the color.
    - Click anywhere on the color spectrum to create a new color bar.

    We recommend setting the color for missing data (-1) to white or light gray
    and using a continuous (not discrete) color scale for the mutation data.

.. _documentation for Fold: https://rna.urmc.rochester.edu/Text/Fold.html
.. _folding constraints file: https://rna.urmc.rochester.edu/Text/File_Formats.html#Constraint
.. _VARNA: https://varna.lisn.upsaclay.fr/
.. _RNAstructure: https://rna.urmc.rochester.edu/RNAstructure.html
.. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/
