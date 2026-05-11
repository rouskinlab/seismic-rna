
Collate: Combine HTML graphs and SVG drawings into a single report
--------------------------------------------------------------------------------

The Collate step gathers HTML graphs (from ``seismic graph``) and SVG structure
drawings (from ``seismic draw``) and combines them into a single, self-contained
HTML report file.
This makes it easy to share or browse all results from a project in one place.

Collate: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass the output directory (or individual HTML/SVG files) as inputs.
SEISMIC-RNA will search recursively for all graphs and drawings::

    seismic collate {out}

Use ``--include-graph`` / ``--no-include-graph`` to control whether HTML graphs
from ``seismic graph`` are included (default: yes).
Use ``--include-svg`` / ``--no-include-svg`` to control whether SVG drawings
from ``seismic draw`` are included (default: yes).

Collate: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Collate setting: Name the report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Use ``--name`` to set the base name for the output HTML file (default:
``collated``).

Use ``--verbose-name`` to append the sample name, reference, region, and graph
type to the file name, which is useful when collating multiple reports in the
same directory.

Collate setting: Group graphs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Use ``--group`` to control how graphs are organized within the report:

=========== ====================================================================
``--group`` Grouping
=========== ====================================================================
``sample``  One section per sample (default)
``graph``   One section per graph type
``branches``One section per analysis branch
``region``  One section per reference region
``all``     All graphs in a single section
=========== ====================================================================

Collate setting: Portable output
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, the output HTML file references the original graph files by their
paths on disk (linked mode).
This means the graphs update automatically if the source files change, but the
HTML file cannot be moved to another computer without the graph files.

Use ``--portable`` to embed all graphs directly into the HTML file.
This makes the file fully self-contained and easy to share, but the file will
be larger and will not reflect future changes to the source graphs.

Collate setting: Output directory
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

By default, the collated HTML file is written to the lowest-level directory
common to all input graphs.
Use ``--collate-out-dir`` to write it to a specific directory instead.

Collate: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Collate step writes a single HTML file, ``{name}.html`` (or a name derived
from ``--name`` and optionally ``--verbose-name``), to the output directory.
Open this file in any web browser to browse all graphs and drawings interactively.
