
General parameters for graphing
--------------------------------------------------------------------------------

Many subcommands of ``seismic graph`` share the same arguments and options.

.. _graph_inputs:

General parameters for graphing: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All graph subcommands accept any number of table files.
See :doc:`../inputs` for ways to list multiple files.

Some subcommands accept only positional tables (``per-pos``), some only per-read
tables (``per-read``), and some both.
The documentation for each subcommand states which table files it accepts.
If you give any table files that are not accepted, then SEISMIC-RNA will simply
ignore them.

.. _graph_data:

General parameters for graphing: Input data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--rels`` (``-r``): Relationship(s) to graph
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For tables where each column is a type of relationship (i.e. ``per-pos`` and
``per-read`` tables), this option selects relationships.
Each type of relationship has a one-letter code:

-   ``v``: Covered
-   ``n``: Informed
-   ``r``: Matched
-   ``m``: Mutated
-   ``s``: Subbed
-   ``a``: Subbed to A
-   ``c``: Subbed to C
-   ``g``: Subbed to G
-   ``t``: Subbed to T
-   ``d``: Deleted
-   ``i``: Inserted

You can give ``-r`` any number of times; each will make its own graph.
For example, ``-r s`` would make one graph:

1. substitutions

And ``-r r -r s -r d`` would make three graphs:

1.  matches
2.  subsititions
3.  deletions

You can also give more than one code after ``-r`` for ``profile`` graphs.
For example, ``-r sdi`` would make one stacked profile graph:

1.  substitutions, deletions, and insertions

And ``-r sdi -r acgt`` would make two stacked profile graphs:

1.  substitutions, deletions, and insertions
2.  subsitutions to A, C, G, and T

``--use-ratio/--use-count``: Graph ratios or counts
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

All tables contain counts.
You can graph either the counts themselves (``--use-count``) or ratios of one
count to another (``--use-ratio``).

-   Graphing counts is useful when you care about the number of reads.
    For example, you can graph the read coverage with ``-r v --use-count``.
-   Graphing ratios is useful when you care about the fraction of reads with a
    certain attribute.
    For example, you can graph mutation rates with ``-r m``.

For ratios, the numerator is the count of the relationship; if that relationship
is Covered or Informed, then the denominator is the count of Covered, otherwise
the count of Informed relationships.
For example, ``-r m`` would use the ratio of Mutated (``m``) to Informed; while
``-r n`` would use the ratio of Informed (``n``) to Covered.

``--quantile`` (``-q``): Quantile for normalizing ratios
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Some graphs (e.g. ``delprof``) work best when the data are normalized.
Use ``-q`` to normalize ratios to a quantile.
See :doc:`../normalize` for more information.
With ``-q``, only ratios are normalized -- not counts.

.. _graph_outputs:

General parameters for graphing: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{sample}/graph/{ref}/{sect}``,
where ``{out}`` is the output directory, ``{sample}`` is the sample, ``{ref}``
is the reference, and ``{sect}`` is the section.
Each output file is prefixed with the name of the subcommand that produced it.

``--cgroup``: Group orders and clusters into output files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When graphing tables containing multiple clusters, you can group

-   ``--cgroup indiv``: Each cluster individually in its own file.
    Each output file will be named with its order and cluster, e.g. ``2-1`` for
    cluster 1 of order 2.
-   ``--cgroup order``: Each order in its own file and every cluster in that
    order as a subplot in that file.
    Each output file will be named with its order, with an "x" for the cluster,
    e.g. ``2-x`` for all clusters of order 2.
-   ``--cgroup unite``: All orders and clusters as subplots in one file.
    Each output file will be named with an "x" for both the order and cluster,
    e.g. ``x-x`` for all clusters of all orders.

``--out-dir`` (``-o``): Destination for all new output files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

All graph files will be written into this output directory.
Graph subcommands need you to specify an output directory because (unlike with
most other commands) the data can come from multiple sources (e.g. table files
from several different output directories, or exogenous RNA structure files).
As a result, where the output files should go could be ambiguous in some cases
if you did not specify an output directory.

Types of graph output files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Each graph subcommand can output graphs in several formats:

-   ``--csv/--no-csv``: **Output the raw data as a CSV file.**
    Useful when you need the graph's raw data, e.g. to analyze further or submit
    along with a paper or figure you are publishing.
-   ``--html/--no-html``: **Output the graph as an interactive HTML file.**
    Useful when you need an interactive graph that you can resize, filter, and
    hover over to see more details.
-   ``--pdf/--no-pdf``: **Output the graph as a PDF file.**
    Useful when you need to edit the graph in a vector graphics program like
    Inkscape or Illustrator.
