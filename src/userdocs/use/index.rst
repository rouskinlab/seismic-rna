********************************************************************************
How to Use
********************************************************************************

This section is the reference manual for every SEISMIC-RNA subcommand.
Each command page describes that command's inputs and outputs, every option,
caveats, performance tips, common errors, and common unexpected results.

Workflow commands are listed in the order in which they are typically run.
Utility commands are independent helpers.
The shared-topic pages at the bottom describe conventions that apply across
multiple commands.


.. toctree::
    :maxdepth: 1
    :caption: Workflow commands

    workflow/wf
    workflow/demult
    workflow/align
    workflow/idmut
    workflow/pool
    workflow/filter
    workflow/cluster
    workflow/join
    workflow/ensembles
    workflow/table
    workflow/fold
    workflow/graph
    workflow/draw
    workflow/collate
    workflow/export


.. toctree::
    :maxdepth: 1
    :caption: Utility commands

    utility/splitbam
    utility/cleanfa
    utility/lists
    utility/renumct
    utility/sim
    utility/importmm
    utility/ct2db
    utility/db2ct
    utility/migrate
    utility/test


.. toctree::
    :maxdepth: 1
    :caption: Shared topics

    inputs
    regions
    logging
    parallel
    branch
