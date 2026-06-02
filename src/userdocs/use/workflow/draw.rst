********************************************************************************
seismic draw
********************************************************************************


Purpose
================================================================================

``seismic draw`` renders predicted RNA secondary structures as 2D diagrams,
coloring each base by its mutation rate reactivity.
It uses RNArtistCore_ to generate publication-quality SVG images.

Requires RNArtistCore_ (installed automatically on first run unless
``--no-update-rnartistcore`` is set).


Inputs
================================================================================

Fold CT files or fold output directories
    One or more CT files or ``fold-report.json`` files produced by
    ``seismic fold``.
    See :doc:`/use/inputs`.

Table CSV files (optional)
    Per-position mutation rate tables to color the bases.
    If omitted, bases are drawn without color.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/draw/{ref}/{reg}/``.

``{profile}-{n}.svg``
    SVG diagram of structure ``n`` of a profile, e.g.
    ``{reg}__average-0.svg`` (one file per drawn structure).

``{profile}-{n}.png``
    PNG diagram (only if ``--draw-png`` is set).


Quick example
================================================================================

Draw all predicted structures for ``sample-1``::

    seismic draw out/sample-1/fold/ref-1


Options
================================================================================

Structure selection
    ``--struct-num N``
        Draw structure number N (0-indexed); ``-1`` for all structures.
        Default: draw the structure with the best AUC-ROC.

Appearance
    ``--color/--no-color``
        Color bases by their reactivity (default on).
    ``--draw-svg/--no-draw-svg``
        Output SVG files (default on).
    ``--draw-png/--no-draw-png``
        Output PNG files (default off).

Maintenance
    ``--update-rnartistcore/--no-update-rnartistcore``
        Check for and install updates to RNArtistCore (default off).

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


See also
================================================================================

- :doc:`fold` — produces the CT files this step draws
- :doc:`table` — produces the reactivity tables used to color bases
- :doc:`/formats/data/ct`
- :doc:`/use/parallel`


.. _RNArtistCore: https://github.com/fjossinet/RNArtistCore
