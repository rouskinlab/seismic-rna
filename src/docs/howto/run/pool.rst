
Pool: Combine samples from the Relate step
--------------------------------------------------------------------------------

Pool: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pool input file: Relate report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of Relate report files as inputs for the Pool step.
See :doc:`../inputs` for ways to list multiple files.

Relate reports will be pooled only if they share both

- the reference
- the top-level output directory, i.e. ``--out-dir`` (``-o``)

For each pair of these two attributes, SEISMIC-RNA will produce a pooled sample
from all Relate reports with those attributes.
The original Relate reports will not be deleted or modified; you will merely get
a new Pool report file for each pooled sample.

For example, if you ran the command ::

    seismic pool -P {pooled} {out}/sample-1/relate/ref-1 {out}/sample-1/relate/ref-2 {out}/sample-2/relate/ref-1 {out}/sample-1/relate/ref-2

where ``{out}`` is the path of your output directory from the Relate step and
``{pooled}`` is the name you want to give to each pooled sample, then you would
get two new Pool reports representing the pooled samples:

- ``{out}/{pooled}/relate/ref-1/relate-report.json``: made from
  ``{out}/sample-1/relate/ref-1/relate-report.json`` and
  ``{out}/sample-2/relate/ref-1/relate-report.json``
- ``{out}/{pooled}/relate/ref-2/relate-report.json``: made from
  ``{out}/sample-1/relate/ref-2/relate-report.json`` and
  ``{out}/sample-2/relate/ref-2/relate-report.json``

To pool all valid combinations of Relate reports in ``{out}`` into samples named
``{pooled}``, you can use the command::

    seismic pool -P {pooled} {out}

Pool: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pool setting: Name of pooled sample
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can choose the name of your pooled sample(s) using ``--pool`` (``-P``).
If you omit this option in ``seismic pool``, then it will default to ``pooled``.
If you omit this option in ``seismic wf``, then the Pool step will not run.

Pool: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{pooled}/relate/{ref}``, where
``{out}`` is the output directory, ``{pooled}`` is the pooled sample name, and
``{ref}`` is the name of the reference.

Pool output file: Pool report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

SEISMIC-RNA writes a Pool report file, ``relate-report.json``, that records the
names of the samples you pooled.
See :doc:`../../formats/report/pool` for more information.
The file is named ``relate-report.json`` not because its contents are identical
to those of a Relate report file (they aren't) but because SEISMIC-RNA can more
easily use Relate and Pool report files interchangably when they have the same
file names.
You can pass both Relate and Pool report files into the Mask and Table steps.

Pool: Troubleshoot and optimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RelateDataset: Field 'Number of Reads' has no default
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This error means that you provided a Pool report file instead of a Relate report
file as input to the Pool step.
You may ignore this error because it will not affect your intended output files.
