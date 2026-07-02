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
- :doc:`workflow/filterscan`
- :doc:`workflow/clusterscan`
- :doc:`utility/splitbam`
- :doc:`utility/importmm`
- :doc:`utility/lists`

``seismic wf`` does not accept ``--branch`` directly, because it runs many steps
at once.  Instead, to branch one or more steps within the workflow, use
``--wf-branch STEP NAME``, giving the name of the step followed by the branch
name.  Repeat the option to branch several steps in one run, e.g.
``--wf-branch filter strict --wf-branch cluster strict``.  The step name must be
one of ``demult``, ``align``, ``idmut``, ``filter``, ``filterscan``, ``cluster``,
``clusterscan``, or ``fold``; any other name raises an error.  See
:doc:`workflow/wf`.


See also
--------------------------------------------------------------------------------

- :doc:`parallel` — parallelize across multiple inputs
- :doc:`regions` — run the same step over different regions of a reference
