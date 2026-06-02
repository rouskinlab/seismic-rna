Parallelize Tasks
================================================================================

SEISMIC-RNA can use multiple CPU cores to speed up processing.


``--num-cpus N``
--------------------------------------------------------------------------------

Pass ``--num-cpus N`` to any command to use up to N cores simultaneously.
By default, SEISMIC-RNA uses all available cores on your machine.

SEISMIC-RNA parallelizes across independent inputs: each core processes a
different sample, reference, or batch at the same time.
The speedup is roughly proportional to N, up to the number of independent
inputs.
For example, if you have 4 references and set ``--num-cpus 4``, all four
are processed simultaneously.

On a shared computing cluster, set ``--num-cpus`` to your allocated core count
rather than leaving it at the default.


.. _batches:

Batches
--------------------------------------------------------------------------------

The IDmut step divides each dataset into batches of reads, controlled by
``--batch-size N`` (number of reads per batch, default 65,536).
Downstream steps (Filter, Cluster) reuse those same batches, so
``--batch-size`` is set only on IDmut (or on :doc:`workflow/wf`, which passes
it through).

Batches serve two purposes:

- **Speed**: multiple batches can be processed simultaneously.
- **Memory**: only one batch needs to fit in RAM at a time.
  Reducing ``--batch-size`` lowers peak memory usage at the cost of more output
  files and slightly more overhead.


Recommendations
--------------------------------------------------------------------------------

- Start with the defaults (all CPUs, batch size 65,536).
- If a step runs out of memory, reduce ``--num-cpus`` or ``--batch-size``.
- If you have many small inputs, set ``--num-cpus`` to match the number of
  inputs for maximum utilization.


See also
--------------------------------------------------------------------------------

- :doc:`branch` — run multiple analyses in parallel with different settings
