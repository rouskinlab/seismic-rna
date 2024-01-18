
Pool: Merge samples (vertically) from the Relate step
--------------------------------------------------------------------------------

Pool: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pool input file: Relate/Pool report
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of Relate report files as inputs for the Pool step.
You can also give Pool report files, to pool samples that were themselves made
by pooling other samples.
See :doc:`../inputs` for ways to list multiple files.

.. note::
    SEISMIC-RNA will not double-count any of your original samples, even if they
    appear in multiple report files you are pooling.
    It will just log a warning if it finds any samples given multiple times.

Relate and Pool reports will be pooled only if they share both

- the top-level output directory, i.e. ``--out-dir`` (``-o``)
- the reference

For each pair of these two attributes, SEISMIC-RNA will produce a pooled sample
from all Relate/Pool reports with those attributes.
The original Relate/Pool report files will not be deleted or modified; you will
merely get a new Pool report file for each pooled sample.

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

To pool all valid combinations of Relate/Pool reports in ``{out}`` into samples
named ``{pooled}``, you can use the command::

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

Pool ... got duplicate samples
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This warning means that an original (unpooled) sample appeared more than once in
the report files you are pooling.

For example, suppose that you pool ``sample-A`` and ``sample-B``::

    seismic pool -P pool-1 out/sample-A out/sample-B

Then you try to pool ``sample-A`` with the pooled sample ``pool-1``::

    seismic pool -P pool-2 out/sample-A out/pool-1

This second command will warn that ``sample-A`` is duplicated because it appears
in both the report files for ``sample-A`` and ``pool-1``.

If you get this warning, then you should check your Pool report file to ensure
it contains all the samples you want and none that you don't.

Overwriting ... would cause data loss
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This error means that you attempted to create a pooled sample with the same name
as an existing non-pooled sample while using ``--force``, e.g. ::

    seismic pool --force -P sample-A out

if ``out/sample-A`` already exists.

Doing so would overwrite the Relate report for the original, non-pooled sample,
making the sample unusable (unless you reran ``seismic relate`` on that sample).
To prevent data loss, the Pool step refuses to overwrite Relate reports, even
with ``--force``.
