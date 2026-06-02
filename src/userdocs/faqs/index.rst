********************************************************************************
FAQs
********************************************************************************

This page collects answers to common questions that come up while installing,
running, and interpreting SEISMIC-RNA.


Installation and external tools
================================================================================

Why does ``seismic fold`` fail with a ``DATAPATH`` error?
--------------------------------------------------------------------------------

The RNAstructure backend needs the ``DATAPATH`` environment variable to
point at its ``data_tables/`` directory (parameter files, free-energy
tables, etc.).
SEISMIC-RNA can usually guess the correct value of ``DATAPATH``, without need
for setting ``DATAPATH`` manually, if RNAstructure was installed using Conda or
Mamba or was downloaded from the RNAstructure website.
Run ``seismic datapath`` to see what SEISMIC-RNA guesses.
If it does not guess ``DATAPATH`` correctly, then set it manually before running
``seismic fold`` or ``seismic wf --fold``::

    export DATAPATH=/whatever/your/path/is/to/RNAstructure/data_tables

Add the line to your shell rc file (``~/.bashrc``, ``~/.zshrc``) to make
it persistent.
See :doc:`/install/index` for the full RNAstructure setup, and
:doc:`/use/workflow/fold` for backend selection.


Why does ``seismic draw`` print "Skip draw if RNARTISTCORE is not installed"?
--------------------------------------------------------------------------------

``seismic draw`` shells out to ``RNArtistCore`` to render structures as
SVG.
If the ``RNARTISTCORE`` environment variable does not point to a working
install, the draw step is skipped with a warning rather than failing the
whole pipeline.
Install RNArtistCore and set::

    export RNARTISTCORE=/path/to/rnartistcore

See :doc:`/install/index` and :doc:`/use/workflow/draw`.


Which optional dependencies do I actually need?
--------------------------------------------------------------------------------

- **Bowtie 2** and **samtools** — required for :doc:`/use/workflow/align`.
- **samtools** — required for :doc:`/use/workflow/idmut`.
- **RNAstructure** (with ``DATAPATH``) — required for the default
  ``--fold-backend rnastructure`` and for ``--pseudoknots``
  (``ShapeKnots``).
- **ViennaRNA** — required for ``--fold-backend viennarna``.
- **RNArtistCore** (with ``RNARTISTCORE``) — required for
  :doc:`/use/workflow/draw`.

If a backend / tool is not installed, the corresponding step fails with a clear
message — earlier steps still work.


Renamed commands and options (v0.25.0)
================================================================================

In v0.25.0 several commands and options were renamed for clarity.
A summary is in `CHANGELOG.md
<https://github.com/rouskinlab/seismic-rna/blob/main/CHANGELOG.md>`_;
:doc:`/use/utility/migrate` migrates older output directories
automatically.

Commands
    +-------------------------+-------------------------+
    | Old name                | New name                |
    +=========================+=========================+
    | ``seismic relate``      | ``seismic idmut``       |
    +-------------------------+-------------------------+
    | ``seismic mask``        | ``seismic filter``      |
    +-------------------------+-------------------------+

Options
    +---------------------------------+--------------------------------------+
    | Old option                      | New option                           |
    +=================================+======================================+
    | ``--fold-table``                | ``--fold-table-region``              |
    +---------------------------------+--------------------------------------+
    | ``--fold-full``                 | ``--fold-full-region``               |
    +---------------------------------+--------------------------------------+
    | ``--quantile``                  | ``--fold-quantile`` (fold) +         |
    |                                 | ``--graph-quantile`` (graph)         |
    +---------------------------------+--------------------------------------+
    | ``--mismatch-tolerence``        | ``--mismatch-tolerance``             |
    +---------------------------------+--------------------------------------+
    | ``--index-tolerence``           | ``--index-tolerance``                |
    +---------------------------------+--------------------------------------+
    | ``--mask-*`` (reads)            | ``--drop-*``                         |
    +---------------------------------+--------------------------------------+
    | ``--mask-gu``                   | ``--mask-a/c/g/u`` (per-base)        |
    +---------------------------------+--------------------------------------+
    | ``--pooled`` (``seismic pool``) | positional pool-name argument        |
    +---------------------------------+--------------------------------------+
    | ``--joined`` (``seismic join``) | positional region-name argument      |
    +---------------------------------+--------------------------------------+

Run ``seismic migrate out -o out_new`` to migrate older output directory
trees to the v0.25.0 layout.
See :doc:`/use/utility/migrate`.


Probes and probe-specific defaults
================================================================================

What does ``--probe`` do?
--------------------------------------------------------------------------------

``--probe {DMS|SHAPE|ETC|none}`` declares which chemical probe was used,
and switches on probe-appropriate defaults across the pipeline:

- **Filter** masks the bases that are unreactive to that probe by default
  (DMS: G + U; ETC: A + C; SHAPE / none: no mask).
- **Filter** picks ``--mut-collisions`` (``drop`` for DMS;
  ``merge`` for SHAPE/ETC).
- **Fold** picks ``--fold-backend`` (RNAstructure for DMS, ViennaRNA for
  SHAPE/ETC/none) and ``--fold-energy-method`` (``cordero`` for DMS,
  ``eddy`` for SHAPE/ETC/none).
- **Sim** turns on ``--min-mut-gap-weights`` for DMS and
  ``--injected-mut-probs`` for SHAPE/ETC by default.

Override any of these explicitly to use a non-default combination.

When should I set ``--probe none``?
--------------------------------------------------------------------------------

For untreated control samples — that is, sequencing reads from RNA that
was not exposed to a chemical probe.
This disables base-specific masking and selects the appropriate folding
defaults.


Performance and parallelism
================================================================================

How do I make ``seismic wf`` run faster?
--------------------------------------------------------------------------------

- **Parallelize.** Use ``--num-cpus N`` (where ``N`` is the number of
  CPUs to use). See :doc:`/use/parallel`.
- **Tune ``--batch-size``.** Smaller batches use less memory per worker
  but produce more files. See :doc:`/use/workflow/idmut`.
- **Disable indel ambiguity.** If your data has few indels, pass
  ``--no-ambindel`` to skip the ambiguous-indel enumeration. See
  :doc:`/use/workflow/idmut`.

How do I limit memory usage?
--------------------------------------------------------------------------------

Reduce ``--batch-size`` and/or ``--num-cpus``.
Monitor memory with ``top`` (Linux/macOS) or Activity Monitor (macOS) during a
run; if memory pressure is high, halve ``--batch-size`` and try again.


Output and interpretation
================================================================================

Where do my output files go?
--------------------------------------------------------------------------------

Under ``--out-dir`` (default ``./out``), organized as
``{out}/{sample}/{step}/{ref}/...``.
Each step also writes a JSON report file (e.g. ``idmut-report.json``,
``filter-report.json``) that records the settings used and summarizes
the results.
See :doc:`/formats/report/index`.


Why is my output empty / why were all my reads filtered out?
--------------------------------------------------------------------------------

The Filter step drops reads that fail any of its criteria.
Open the filter report (``filter-report.json``) to see how many reads
each criterion dropped.
Common causes:

- ``--min-fcov-read`` set too high (drops reads that do not cover the
  full region).
- ``--probe`` set incorrectly (wrong default base masking).
- ``--mut-collisions drop`` dropping reads with closely-spaced
  mutations.
- The IDmut step itself produced empty batches because the alignment
  maps had too few reads (see :doc:`/use/workflow/idmut`).


I ran ``seismic wf`` but got no structures / no clustered output / no JSON.
--------------------------------------------------------------------------------

The Fold, Cluster, Draw, and Export steps are **off by default** in
``seismic wf``.
Turn them on explicitly:

- Cluster: ``-k N`` (or ``--max-clusters N``)
- Fold: ``--fold``
- Draw: ``--draw``
- Export: ``--export``

See :doc:`/use/workflow/wf`.


How do I re-run only one step?
--------------------------------------------------------------------------------

Run the step's command directly on the upstream outputs, e.g. ::

    seismic filter out/sampleA/idmut/refX/idmut-report.json --probe DMS

``seismic wf`` will also skip any step whose outputs already exist;
pass ``--force`` to overwrite.


Troubleshooting and support
================================================================================

How do I see what SEISMIC-RNA is actually doing?
--------------------------------------------------------------------------------

Use ``-v`` / ``-vv`` / ``-vvv`` / ``-vvvv`` to increase log verbosity,
and ``-q`` / ``-qq`` / ``-qqq`` / ``-qqqq`` to decrease it.
All messages are also written to a log file regardless of verbosity.
See :doc:`/use/logging`.

I found a bug. Where do I report it?
--------------------------------------------------------------------------------

See :doc:`/issues`.

I have a question that's not answered here.
--------------------------------------------------------------------------------

Open a GitHub issue (see :doc:`/issues`) or check the
:doc:`/tutorials/index` for worked examples.
