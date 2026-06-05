
Log Messages
================================================================================

As SEISMIC-RNA runs, it writes messages to the terminal (standard error) and to
a log file.
Messages are ranked by level from most to least severe:

======= ========== =============================================================
 Level   Name       Used for
======= ========== =============================================================
 -2      ERROR      Problem preventing one output
 -1      WARNING    Abnormal event that is fully recovered
  0      INFO       Current pipeline step (sample / reference / region)
  1      DEBUG      Repeated subtasks: batch I/O, directory scans, shell commands
  2      TRACE      Fine-grained detail: internal variable values
======= ========== =============================================================


Controlling terminal output
--------------------------------------------------------------------------------

By default, INFO (0) and more-severe messages appear in the terminal.
Use ``--verbose`` (``-v``) to see more; ``--quiet`` (``-q``) to see less.
Each additional ``v`` or ``q`` shifts the threshold by one level:

====== ===========================================================================
 Flag   Shows in terminal
====== ===========================================================================
 -vv    TRACE and above (everything)
 -v     DEBUG and above
 (none) INFO and above (default)
 -q     WARNING and above
 -qq    ERROR only
====== ===========================================================================

These flags go between ``seismic`` and the subcommand::

    seismic -v cluster -k 3 out/sars2-fse


Log files
--------------------------------------------------------------------------------

By default, messages are also written to a log file at::

    ./log/seismic-rna_YYYY-MM-DD_hh-mm-ss.log

Use ``--log PATH`` (between ``seismic`` and the subcommand) to choose a
different path, or ``--log ""`` to disable log files.

The log file records messages at **the same verbosity level as the terminal**,
so ``-v`` enlarges both the terminal output and the log file.
Each line in the log file is prefixed with a timestamp::

    2026-06-04 14:22:01 INFO    12345  Clustering sars2-fse / region fse

To watch a log file update in real time, open it with ``less`` and press
Shift-F::

    less log/seismic-rna_2024-04-08_15-20-09.log


.. _standard error: https://en.wikipedia.org/wiki/Standard_streams#Standard_error_(stderr)
.. _less: https://greenwoodsoftware.com/less/
