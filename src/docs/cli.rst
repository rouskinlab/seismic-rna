********************************************************************************
Commands, Arguments, Options
********************************************************************************

.. note::
    All subcommands, arguments, and options listed in these documents can be
    printed, along with a short explanation of each, on the command line by
    running ``seismic [command(s)] --help``.
    For example,
    ``seismic --help`` prints all options and subcommands for the command
    ``seismic``,
    ``seismic graph --help`` prints all options and subcommands for the command
    ``seismic graph``, and
    ``seismic graph profile --help`` prints all options and subcommands for the
    command ``seismic graph profile``.

Run all steps of the workflow
================================================================================

.. _cli_wf:

.. click:: seismicrna.workflow:cli
    :prog: seismic wf

Run individual steps of the workflow
================================================================================

.. _cli_demult:

.. click:: seismicrna.demult:cli
    :prog: seismic demult

.. _cli_align:

.. click:: seismicrna.align:cli
    :prog: seismic align

.. _cli_relate:

.. click:: seismicrna.relate:cli
    :prog: seismic relate

.. _cli_mask:

.. click:: seismicrna.mask:cli
    :prog: seismic mask

.. _cli_cluster:

.. click:: seismicrna.clust:cli
    :prog: seismic cluster

.. _cli_table:

.. click:: seismicrna.table:cli
    :prog: seismic table

.. _cli_fold:

.. click:: seismicrna.fold:cli
    :prog: seismic fold

.. _cli_graph:

seismic graph
--------------------------------------------------------------------------------

.. click:: seismicrna.graph.profile:cli
    :prog: seismic graph profile

.. click:: seismicrna.graph.delprof:cli
    :prog: seismic graph delprof

.. click:: seismicrna.graph.corroll:cli
    :prog: seismic graph corroll

.. click:: seismicrna.graph.scatter:cli
    :prog: seismic graph scatter

Extra Utilities
================================================================================

.. note::
    For every extra utility (that is not part of the main workflow), the name
    begins with ``+``.

.. _cli_cleanfa:

.. click:: seismicrna.cleanfa:cli
    :prog: seismic +cleanfa

.. _cli_renumct:

.. click:: seismicrna.renumct:cli
    :prog: seismic +renumct

.. _cli_test:

.. click:: seismicrna.test:cli
    :prog: seismic +test