********************************************************************************
Command Line Reference
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

.. _cli_idmut:

.. click:: seismicrna.idmut:cli
    :prog: seismic idmut

.. _cli_pool:

.. click:: seismicrna.pool:cli
    :prog: seismic pool

.. _cli_filter:

.. click:: seismicrna.filter:cli
    :prog: seismic filter

.. _cli_cluster:

.. click:: seismicrna.cluster:cli
    :prog: seismic cluster

.. _cli_join:

.. click:: seismicrna.join:cli
    :prog: seismic join

.. _cli_ensembles:

.. click:: seismicrna.ensembles.main:cli
    :prog: seismic ensembles

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

.. click:: seismicrna.graph.giniroll:cli
    :prog: seismic graph giniroll

.. click:: seismicrna.graph.snrroll:cli
    :prog: seismic graph snrroll

.. click:: seismicrna.graph.histpos:cli
    :prog: seismic graph histpos

.. click:: seismicrna.graph.histread:cli
    :prog: seismic graph histread

.. click:: seismicrna.graph.mutdist:cli
    :prog: seismic graph mutdist

.. click:: seismicrna.graph.poscorr:cli
    :prog: seismic graph poscorr

.. click:: seismicrna.graph.roc:cli
    :prog: seismic graph roc

.. click:: seismicrna.graph.aucroll:cli
    :prog: seismic graph aucroll

.. click:: seismicrna.graph.abundance:cli
    :prog: seismic graph abundance

.. _cli_draw:

.. click:: seismicrna.draw.main:cli
    :prog: seismic draw

.. _cli_collate:

.. click:: seismicrna.collate:cli
    :prog: seismic collate

.. _cli_export:

.. click:: seismicrna.export.main:cli
    :prog: seismic export

Utility commands
================================================================================

.. _cli_splitbam:

.. click:: seismicrna.align.split:cli
    :prog: seismic splitbam

.. _cli_cleanfa:

.. click:: seismicrna.cleanfa:cli
    :prog: seismic cleanfa

.. _cli_list:

.. click:: seismicrna.lists:cli
    :prog: seismic list

.. _cli_renumct:

.. click:: seismicrna.renumct:cli
    :prog: seismic renumct

.. _cli_ct2db:

.. click:: seismicrna.core.rna.convert:cli_ct2db
    :prog: seismic ct2db

.. _cli_db2ct:

.. click:: seismicrna.core.rna.convert:cli_db2ct
    :prog: seismic db2ct

.. _cli_importmm:

.. click:: seismicrna.importmm.main:cli
    :prog: seismic importmm

.. _cli_migrate:

.. click:: seismicrna.migrate:cli
    :prog: seismic migrate

.. _cli_fold_datapath:

.. click:: seismicrna.fold.datapath:cli_datapath
    :prog: seismic fold datapath

.. _cli_test:

.. click:: seismicrna.test:cli
    :prog: seismic test

seismic sim
--------------------------------------------------------------------------------

.. click:: seismicrna.sim.total:cli
    :prog: seismic sim total

.. click:: seismicrna.sim.ref:cli
    :prog: seismic sim ref

.. click:: seismicrna.sim.fold:cli
    :prog: seismic sim fold

.. click:: seismicrna.sim.params:cli
    :prog: seismic sim params

.. click:: seismicrna.sim.muts:cli
    :prog: seismic sim muts

.. click:: seismicrna.sim.ends:cli
    :prog: seismic sim ends

.. click:: seismicrna.sim.clusts:cli
    :prog: seismic sim clusts

.. click:: seismicrna.sim.idmut:cli
    :prog: seismic sim idmut

.. click:: seismicrna.sim.fastq:cli
    :prog: seismic sim fastq

.. click:: seismicrna.sim.abstract:cli
    :prog: seismic sim abstract
