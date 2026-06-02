
General parameters for graphing
--------------------------------------------------------------------------------

Many subcommands of ``seismic graph`` share the same arguments and options.

.. _graph_inputs:

General parameters for graphing: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All graph subcommands accept any number of table files.
See :doc:`/use/inputs` for ways to list multiple files.

Some subcommands accept only per-position tables, some only per-read
tables, and some both.
The documentation for each subcommand states which table files it accepts.
If you give any table files that are not accepted, then SEISMIC-RNA will simply
ignore them.

.. _graph_data:

General parameters for graphing: Input data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--rels`` (``-r``): Relationship(s) to graph
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For tables where each column is a type of relationship (i.e. per-position and
per-read tables), this option selects relationships.
Each type of relationship has a one-letter code:

-   ``v``: Covered
-   ``n``: Informative
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

Using a ratio divides the count of the relationship by the count of Covered if
that relationship is Covered or Informative, otherwise Informative.
For example, ``-r m`` would use the ratio of Mutated (``m``) to Informative;
while ``-r n`` would use the ratio of Informative (``n``) to Covered.

``--graph-quantile``: Quantile for normalizing ratios
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Some graphs (e.g. ``delprof``) work best when the data are normalized.
Use ``--graph-quantile`` to normalize ratios to a quantile.
See :doc:`/algos/normalize` for more information.
With ``--graph-quantile``, only ratios are normalized -- not counts.

.. _graph_outputs:

General parameters for graphing: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{sample}/graph/{ref}/{reg}``,
where ``{out}`` is the output directory, ``{sample}`` is the sample, ``{ref}``
is the reference, and ``{reg}`` is the region.
Each output file is prefixed with the name of the subcommand that produced it.

``--cgroup``: Group clusters into output files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When graphing tables containing multiple clusters, you can group

-   ``--cgroup c``: Each cluster individually in its own file.
    Each output file will be named with its number of clusters (K) and cluster,
    e.g. ``2-1`` for cluster 1 of K=2.
-   ``--cgroup k``: Each K in its own file and every cluster of that K as a
    subplot in that file.
    Each output file will be named with its K, with an "x" for the cluster,
    e.g. ``2-x`` for all clusters of K=2.
-   ``--cgroup a``: All clusters (of all Ks) as subplots in one file.
    Each output file will be named with an "x" for both the K and cluster,
    e.g. ``x-x`` for all clusters of all Ks.

``--out-dir`` (``-o``): Destination for combined output files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Subcommands that combine multiple input sources (such as ``scatter`` and
``corroll``, which compare two profiles) accept ``--out-dir`` to say where the
combined output should go, because the data can come from several different
output directories (or exogenous RNA structure files), which would otherwise
make the destination ambiguous.
Single-source subcommands (such as ``profile``) instead write next to their
input data, under ``{out}/{sample}/graph/{ref}/{reg}``.

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
