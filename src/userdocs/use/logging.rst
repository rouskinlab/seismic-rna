
Log Messages
================================================================================

As SEISMIC-RNA runs, it writes messages to a log file and to the terminal
(standard error).
Messages are ranked by level from most to least severe:

======= ========== =============================================================
 Level   Name       Used for
======= ========== =============================================================
 -3      FATAL      Problem preventing all outputs
 -2      ERROR      Problem preventing one output
 -1      WARNING    Abnormal event that is fully recovered
  0      STATUS     Current command or workflow step
  1      TASK       Individual parallelizable task
  2      ACTION     File written or shell command run
  3      ROUTINE    Internal function call
  4      DETAIL     Fine-grained detail
======= ========== =============================================================


Controlling terminal output
--------------------------------------------------------------------------------

By default, STATUS (0) and more-severe messages appear in the terminal.
Use ``--verbose`` (``-v``) to see more; ``--quiet`` (``-q``) to see less.
Each additional ``v`` or ``q`` shifts the threshold by one level:

======= =========================================================================
 Flag    Shows in terminal
======= =========================================================================
 -vvvv   DETAIL and above (everything)
 -vvv    ROUTINE and above
 -vv     ACTION and above
 -v      TASK and above
 (none)  STATUS and above (default)
 -q      WARNING and above
 -qq     ERROR and above
 -qqq    FATAL only
 -qqqq   Nothing
======= =========================================================================

These flags go between ``seismic`` and the subcommand::

    seismic -vv cluster -k 3 out/sars2-fse


Log files
--------------------------------------------------------------------------------

Every message (regardless of terminal verbosity) is written to a log file at::

    ./log/seismic-rna_YYYY-MM-DD_hh-mm-ss.log

Use ``--log PATH`` (between ``seismic`` and the subcommand) to choose a
different path, or ``--log ""`` to disable log files.

To watch a log file update in real time, open it with ``less`` and press
Shift-F::

    less log/seismic-rna_2024-04-08_15-20-09.log


.. _standard error: https://en.wikipedia.org/wiki/Standard_streams#Standard_error_(stderr)
.. _less: https://greenwoodsoftware.com/less/
