********************************************************************************
seismic filter
********************************************************************************

.. note::

    ``filter`` was named ``mask`` in versions before v0.25.0.
    See the project ``CHANGELOG.md`` for the full rename list.


Purpose
================================================================================

``seismic filter`` prepares IDmut output for downstream analysis.
It decides which base changes (substitutions, deletions, insertions) count as
mutations, masks positions unsuitable for your chemical probe (e.g. G and U
for DMS), and drops reads that fail quality thresholds such as too few covered
positions or too high a mutation rate.
The result is a cleaned dataset ready for quantification and structure prediction.


Inputs
================================================================================

IDmut output directory or report file
    One or more IDmut output directories or ``idmut-report.json`` files.
    See :doc:`/use/inputs` for ways to select multiple inputs at once.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/filter/{ref}/{reg}/``.

``filter-batch-{num}.brickle``
    One batch of filtered reads as a ``FilterBatchIO`` object.
    See :doc:`/formats/data/brickle`.

``filter-report.json``
    Summary of settings and results (reads/positions kept and dropped per filter).
    See :doc:`/formats/report/filter`.


Quick example
================================================================================

Filter all IDmut output for ``sample-1`` against ``ref-1`` with DMS defaults::

    seismic filter out/sample-1/idmut/ref-1


Options
================================================================================

Probe preset
    ``--probe {DMS|SHAPE|ETC|none}``
        Sets base masking and collision defaults.
        ``DMS`` (default): masks G and T/U, uses ``--mut-collisions drop``.
        ``ETC``: masks A and C, uses ``--mut-collisions merge``.
        ``SHAPE``: keeps all bases, uses ``--mut-collisions merge``.
        ``none``: keeps all bases and sets ``--min-mut-gap 0`` (so no collision
        handling is applied).
        Individual options override the preset.

Regions
    Select a sub-range of each reference with ``--region-coords`` or
    ``--region-primers``.
    See :doc:`/use/regions`.

Mutation definitions
    ``--count-del/--no-del`` / ``--count-ins/--no-ins``
        Count deletions/insertions as mutations (on by default).
    ``--no-mut XY``
        Ignore a specific mutation type (e.g. ``ag`` = A→G substitution).
    ``--only-mut XY``
        Count only this type as a mutation (overrides all other mutation settings).

Mask out positions
    ``--mask-a`` / ``--mask-c`` / ``--mask-g`` / ``--mask-u``
        Mask/keep positions by base identity.
        G and T/U masked for DMS by default; A and C for ETC; all kept otherwise.
    ``--mask-polya N``
        Mask positions inside poly(A) runs of length ≥ N (default 5; ``0`` to disable).
    ``--mask-pos REF POS`` / ``--mask-pos-file FILE``
        Mask specific positions or positions from a :doc:`/formats/list/listpos` file.

Filter reads
    ``--min-ncov-read N``
        Minimum positions a read must cover (default 1).
    ``--min-fcov-read F``
        Minimum fraction of the region a read must cover (default 0.0).
    ``--min-finfo-read F``
        Minimum fraction of informative positions per read (default 0.95).
    ``--max-fmut-read F``
        Maximum fraction of mutated positions per read (default 1.0; disabled).
    ``--drop-discontig/--keep-discontig``
        Drop paired-end reads with non-contiguous mates (default: drop).
    ``--min-mut-gap N``
        Minimum gap between two mutations in a read; pairs closer than N bases
        are handled by ``--mut-collisions``.
        Probe-specific defaults: 4 for DMS, 2 for SHAPE/ETC, 0 for none.
    ``--mut-collisions {auto|drop|merge}``
        How to handle reads that violate ``--min-mut-gap`` (default ``auto``).

Filter positions
    ``--min-ninfo-pos N``
        Minimum informative reads at a position to keep it (default 1000).
    ``--max-fmut-pos F``
        Maximum fraction of mutated reads at a position (default 1.0; disabled).

Branches
    ``--branch X``
        Create a new branch: output results in ``{out}/{sample}/filter_{branch}``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

Too many reads filtered out
    Check ``filter-report.json`` to see which criterion removed the most reads.
    Loosen ``--min-finfo-read``, ``--min-mut-gap``, or ``--max-fmut-read``, or
    shorten the region so reads cover it fully.
    Excessive mutations may indicate low RNA quality, bad base calls, or
    misalignment.

Too many positions filtered out
    Check ``filter-report.json``.
    Loosen ``--min-ninfo-pos`` or ``--max-fmut-pos``.
    Too few informative reads usually means low sequencing depth or too many
    reads were removed first.

Filter hangs or crashes early
    Almost always out of memory.
    Reduce ``--num-cpus`` or rerun IDmut with a smaller ``--batch-size``.


See also
================================================================================

- :doc:`idmut` — produces the output this step consumes
- :doc:`table` — computes per-position mutation rates from filter output
- :doc:`/formats/data/brickle`, :doc:`/formats/report/filter`
- :doc:`/use/inputs`, :doc:`/use/regions`, :doc:`/use/parallel`
