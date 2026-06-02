********************************************************************************
seismic collate
********************************************************************************


Purpose
================================================================================

``seismic collate`` assembles graphs and structure drawings from multiple
samples or regions into a single interactive HTML report.
Use it at the end of your analysis to view all results side by side in a
web browser.


Inputs
================================================================================

Graph HTML files or draw SVG files, or the directories containing them
    Pass the output directories from ``seismic graph`` and/or ``seismic draw``.
    See :doc:`/use/inputs`.


Outputs
================================================================================

``{name}.html``
    Single interactive HTML report containing all specified graphs and drawings.
    Written to the lowest directory common to all inputs (or ``--collate-out-dir``
    if specified).


Quick example
================================================================================

Collate all graphs and drawings for a sample::

    seismic collate out/sample-1/graph out/sample-1/draw


Options
================================================================================

Report name
    ``--name NAME``
        Prefix of the output HTML file (default ``collated``).
    ``--verbose-name/--no-verbose-name``
        Add file path information to the report name (default off).
    ``--collate-out-dir DIR``
        Write the report to this directory instead of the common ancestor
        of the inputs.

Content
    ``--include-graph/--no-include-graph``
        Include graphs from ``seismic graph`` (default on).
    ``--include-svg/--no-include-svg``
        Include structure drawings from ``seismic draw`` (default on).

Layout
    ``--group {sample|graph|branches|region|all}``
        Group the included plots by this attribute (default ``sample``).
    ``--portable/--no-portable``
        Embed all content directly in the HTML for easy sharing (default off).
        Produces a larger file but works without an internet connection.

Other
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`graph` — produces the graphs included in the report
- :doc:`draw` — produces the structure diagrams included in the report
- :doc:`/use/inputs`
