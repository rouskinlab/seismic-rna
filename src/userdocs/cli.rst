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

Run the entire workflow
================================================================================

.. _cli_wf:

.. click:: seismicrna.wf:cli
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

.. click:: seismicrna.cluster:cli
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

.. click:: seismicrna.graph.scatter:cli
    :prog: seismic graph scatter

.. click:: seismicrna.graph.corroll:cli
    :prog: seismic graph corroll

.. click:: seismicrna.graph.histpos:cli
    :prog: seismic graph histpos

.. click:: seismicrna.graph.histread:cli
    :prog: seismic graph histread

.. click:: seismicrna.graph.roc:cli
    :prog: seismic graph roc

.. click:: seismicrna.graph.aucroll:cli
    :prog: seismic graph aucroll

Extra Utilities
================================================================================

.. note::
    For every extra utility (that is not part of the main workflow), the name
    begins with ``+``.

.. _cli_listpos:

.. click:: seismicrna.lists.listpos:cli
    :prog: seismic +listpos

.. _cli_addclust:

.. click:: seismicrna.cluster.addclust:cli
    :prog: seismic +addclust

.. _cli_delclust:

.. click:: seismicrna.cluster.delclust:cli
    :prog: seismic +delclust

.. _cli_cleanfa:

.. click:: seismicrna.cleanfa:cli
    :prog: seismic +cleanfa

.. _cli_renumct:

.. click:: seismicrna.renumct:cli
    :prog: seismic +renumct

.. _cli_test:

.. click:: seismicrna.test:cli
    :prog: seismic +test

seismic +sim
--------------------------------------------------------------------------------

.. click:: seismicrna.sim.total:cli
    :prog: seismic +sim total

.. click:: seismicrna.sim.ref:cli
    :prog: seismic +sim ref

.. click:: seismicrna.sim.fold:cli
    :prog: seismic +sim fold

.. click:: seismicrna.sim.params:cli
    :prog: seismic +sim params

.. click:: seismicrna.sim.muts:cli
    :prog: seismic +sim muts

.. click:: seismicrna.sim.ends:cli
    :prog: seismic +sim ends

.. click:: seismicrna.sim.clusts:cli
    :prog: seismic +sim clusts
   
.. click:: seismicrna.sim.relate:cli
    :prog: seismic +sim relate

.. click:: seismicrna.sim.fastq:cli
    :prog: seismic +sim fastq

