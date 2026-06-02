********************************************************************************
seismic join
********************************************************************************


Purpose
================================================================================

``seismic join`` concatenates two or more regions of a reference into a single
joined region.
Use it when you defined separate regions (e.g. two amplicons on the same RNA)
and want to analyze or visualize them as one continuous stretch.
The original regions are not modified; join creates a new virtual region.


Inputs
================================================================================

Filter, Cluster, or Join output directories or report files
    One or more paths to the regions you want to join.
    All inputs must belong to the same sample and reference.
    See :doc:`/use/inputs`.

Joined region name
    First positional argument: the name to give the new joined region.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/{step}/{ref}/{joined-region}/``,
where ``{step}`` is ``filter`` or ``cluster`` depending on the inputs.

``filter-report.json`` or ``cluster-report.json``
    Join report describing which regions were combined.
    See :doc:`/formats/report/join`.
    Named to match its source so downstream steps (Table) accept it
    interchangeably with Filter/Cluster reports.


Quick example
================================================================================

Join two amplicon regions into one::

    seismic join joined out/sample-1/filter/ref-1/amp-1 out/sample-1/filter/ref-1/amp-2

To join all regions in the output directory::

    seismic join joined out/


Options
================================================================================

Cluster assignment
    ``--join-clusts FILE`` (``-j``)
        CSV file specifying which cluster of each input region maps to which
        cluster of the joined region (optional).
        See :doc:`/formats/meta/joined`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

Warning: duplicate regions
    An original region appears more than once among the inputs.
    This is harmless — it will not be double-counted — but check that the
    resulting joined region contains exactly what you intended.

Error: would overwrite existing non-joined region
    You chose a name that already exists as a real region.
    Choose a different name or remove the conflicting directory first.


See also
================================================================================

- :doc:`filter` / :doc:`cluster` — produce the data this step joins
- :doc:`table` — accepts join output alongside filter/cluster output
- :doc:`/formats/report/join`, :doc:`/formats/meta/joined`
- :doc:`/use/inputs`, :doc:`/use/parallel`
