************************************************************************
Command Line Interface
************************************************************************

These documents list all commands and their parameters.

They are meant to be complete yet concise references for the commands
and their parameters, but not to provide advice on using them.
For guides on processing your own data, see :doc:`../manuals/index`.

They are meant to be browsed freely, not followed from start to finish.
For step-by-step tutorials, see :doc:`../tutorials/index`.


.. note::

  All subcommands, arguments, and options listed in these documents can
  be printed, along with a short explanation of each, on the command
  line by running ``seismic [command(s)] --help``.
  For example,
  ``seismic --help`` prints all options and subcommands for the command
  ``seismic``,
  ``seismic graph --help`` prints all options and subcommands for the
  command ``seismic graph``,
  and ``seismic graph seqbar --help`` prints all options and subcommands
  for the command ``seismic graph seqbar``.


.. _cli_all:

.. click:: seismicrna.workflow:cli
    :prog: seismic all

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
========================================================================

.. click:: seismicrna.graph.seqbar:cli
    :prog: seismic graph seqbar

.. click:: seismicrna.graph.seqdiff:cli
    :prog: seismic graph seqdiff

.. click:: seismicrna.graph.seqcorr:cli
    :prog: seismic graph seqcorr

.. click:: seismicrna.graph.scatter:cli
    :prog: seismic graph scatter

.. _cli_fastaclean:

.. click:: seismicrna.fastaclean:cli
    :prog: seismic fastaclean
