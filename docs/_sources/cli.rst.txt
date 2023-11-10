.. _cli:

************************************************************************
Command Line Interface
************************************************************************

Target audience: intermediate and experienced users.

These documents detail the arguments and options of each command.
They are references to assist with running SEISMIC-RNA on the command
line, not guides for using each feature of SEISMIC-RNA.

- For guides on using each component of SEISMIC-RNA, see :ref:`manuals`
  (target audience: novice and intermediate users).
- For step-by-step tutorials, see :ref:`tutorials`
  (target audience: novice users).

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

.. click:: seismicrna.main:all_cli
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
========================================================================

.. click:: seismicrna.graph.seqbar:cli
    :prog: seismic graph seqbar

.. click:: seismicrna.graph.seqdiff:cli
    :prog: seismic graph seqdiff

.. click:: seismicrna.graph.seqcorr:cli
    :prog: seismic graph seqcorr

.. click:: seismicrna.graph.scatter:cli
    :prog: seismic graph scatter
