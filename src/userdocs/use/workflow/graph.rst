********************************************************************************
seismic graph
********************************************************************************


Purpose
================================================================================

``seismic graph`` creates plots from mutation rate tables and/or predicted
structures.
It is a group of subcommands, each producing a different plot type.


Subcommands
================================================================================

Run ``seismic graph {subcommand} --help`` for each subcommand's full option list.

============= ===================================================================
Subcommand    Description
============= ===================================================================
``profile``   Bar graph of mutation rate (or other relationship) per position
``delprof``   Bar graph of the difference between two profiles per position
``scatter``   Scatter plot comparing two profiles
``corroll``   Rolling correlation between two profiles
``giniroll``  Rolling Gini coefficient along a profile
``snrroll``   Rolling signal-to-noise ratio along a profile
``aucroll``   Rolling AUC-ROC comparing a profile to a predicted structure
``histpos``   Histogram of relationship values across positions
``histread``  Histogram of relationship values across reads
``mutdist``   Distribution of distances between mutations within reads
``poscorr``   Correlation between pairs of positions
``roc``       ROC curve comparing a mutation rate profile to a structure
``abundance`` Bar graph of cluster abundances
============= ===================================================================


Inputs
================================================================================

Table CSV files, fold CT files, or output directories
    Subcommands that compare two datasets accept two sets of inputs.
    See :doc:`/use/inputs`.


Outputs
================================================================================

``{graph-type}.html``
    Interactive HTML plot (viewable in any web browser).
``{graph-type}.svg`` or ``{graph-type}.png``
    Static image (depending on options).

Outputs go into ``{out}/{sample}/graph/{ref}/{reg}/``.


Quick example
================================================================================

Plot mutation rate profiles for all filter outputs::

    seismic graph profile out/

Compare two profiles::

    seismic graph scatter out/sample-1/filter/ref-1 out/sample-2/filter/ref-1


Options
================================================================================

Common options shared by all graph subcommands include:

``--rels CODES`` (``-r``)
    Which relationship(s) to graph (e.g. ``m`` for mutated, ``n`` for
    unambiguous; default ``m``).
``--use-ratio`` / ``--use-count``
    Graph ratios (e.g. mutation rates) or raw counts (default: ``--use-ratio``).
``--graph-quantile F``
    Normalize and winsorize ratios to this quantile before graphing
    (default 0.0, i.e. no normalization).
``--cgroup {c|k|a}``
    Put each cluster in its own file (``c``), each number of clusters K in its
    own file (``k``), or all clusters in one file (``a``).
``--html`` / ``--no-html``
    Write an interactive HTML plot (default on).
``--svg`` / ``--png`` / ``--pdf``
    Also write static SVG, PNG, and/or PDF images (default off).
``--csv`` / ``--no-csv``
    Also write the graphed data as a CSV file (default on).
``--force``
    Overwrite existing outputs.

See the auto-generated :doc:`/cli` for the complete per-subcommand option lists.


See also
================================================================================

- :doc:`table` — produces the mutation rate tables used as input
- :doc:`fold` — produces the CT files used by ROC/AUC subcommands
- :doc:`/use/inputs`, :doc:`/use/parallel`


.. toctree::
    :hidden:

    graph_assets/index
