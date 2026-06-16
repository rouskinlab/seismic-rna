
Log Messages
================================================================================

As SEISMIC-RNA runs, it writes messages to the terminal (standard error) and to
a log file.
Messages are ranked by level from most to least important:

======= ========== =============================================================
 Level   Name       Used for
======= ========== =============================================================
 -2      Error      Problem causing all or part of SEISMIC-RNA to fail
 -1      Warning    Abnormal event from which recovery is possible
  0      Info       Major normal events: helpful for monitoring progress
  1      Debug      Minor normal events: helpful for monitoring more closely
  2      Trace      Finer detail: helpful for troubleshooting
======= ========== =============================================================


Controlling terminal output
--------------------------------------------------------------------------------

By default, Info (0) and more important messages appear in the terminal.
Use ``--verbose`` (``-v``) to see more; ``--quiet`` (``-q``) to see less.
Each additional ``-v`` or ``-q`` shifts the threshold by one level:

======== =======================================================================
 Flag     Shows in terminal
======== =======================================================================
 -vv      Trace and above (everything)
 -v       Debug and above
 (none)   Info and above (default)
 -q       Warning and above
 -qq      Error only
======== =======================================================================

These flags must go between ``seismic`` and the subcommand::

    seismic -v cluster out/sars2-fse


Log files
--------------------------------------------------------------------------------

By default, messages are also written to a log file at::

    ./log/seismic-rna_YYYY-MM-DD_hh-mm-ss.log

Use ``--log PATH`` (between ``seismic`` and the subcommand) to choose a
different path, or ``--log ""`` to disable log files.

The log file records messages at **the same verbosity level as the terminal**,
so ``-v`` enlarges both the terminal output and the log file.
Each line in the log file is prefixed with a timestamp::

    2026-06-04 14:22:01.123456  Info      12345 Clustering sars2-fse / region fse

To view a log file, open it with less_:

    less -R log/seismic-rna_2024-04-08_15-20-09.log

The ``-R`` flag causes less_ to render the document in color (if the default
colored log messages are enabled); without it, less_ just prints the ANSI codes
for the colors, which makes the document harder to read.

To see new messages appear instantly as they are written to the log file, press
Shift-F when less_ is open.
Press Ctrl-C to stop the continuous updates, then q to quit less_.

.. _standard error: https://en.wikipedia.org/wiki/Standard_streams#Standard_error_(stderr)
.. _less: https://greenwoodsoftware.com/less/
