Branches
================================================================================

The ``--branch`` option (``-b``) lets you run a step under a different name so
its outputs go into a separate directory alongside the default outputs.
This is useful when you want to try different settings on the same data without
overwriting existing results.


How branches work
--------------------------------------------------------------------------------

When you set ``--branch NAME``, the step writes its outputs to
``{out}/{sample}/{step}_{NAME}/`` instead of ``{out}/{sample}/{step}/``.
For example::

    seismic filter --branch strict --min-finfo-read 0.99 out/sample-1/idmut/ref-1

writes to ``out/sample-1/filter_strict/ref-1/``, leaving the default
``out/sample-1/filter/ref-1/`` untouched.

Downstream steps can then run on either branch independently::

    seismic cluster out/sample-1/filter_strict/ref-1
    seismic cluster out/sample-1/filter/ref-1


Which steps support ``--branch``
--------------------------------------------------------------------------------

- :doc:`workflow/demult`
- :doc:`workflow/align`
- :doc:`workflow/idmut`
- :doc:`workflow/filter`
- :doc:`workflow/cluster`
- :doc:`workflow/fold`
- :doc:`workflow/ensembles`
- :doc:`utility/splitbam`
- :doc:`utility/importmm`
- :doc:`utility/lists`

``seismic wf`` does not accept ``--branch``; run individual steps when you need
branched analyses.


See also
--------------------------------------------------------------------------------

- :doc:`parallel` — parallelize across multiple inputs
- :doc:`regions` — run the same step over different regions of a reference
