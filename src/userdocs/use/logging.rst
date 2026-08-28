
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


Progress bars
--------------------------------------------------------------------------------

While a task runs, SEISMIC-RNA shows a progress bar at the bottom of the
terminal, such as::

    align                 |█████         | step 2/8   [00:42<02:31]
    aligning reads        |████████▌     | sample 3/7 [00:42<01:38]

Each bar names the task, what it is counting and which one it is on, how long
it has been running, and roughly how much longer it will take.
The numbers always say what they count, so ``sample 3/7`` can only mean the
3rd of 7 samples.
Everything is numbered from 1, so these bars mean that the workflow is on step
2 of 8 (``align``), which is aligning the 3rd of its 7 samples.
The gauge always matches the numbers beside it: a task on its last item shows
a full gauge, as does a task with only one item.

Tasks nest, and each one gets its own line, the task that contains the others
on top.
The names are padded to the same width so that the gauges line up in a column,
and the whole block stays at the bottom of the terminal.
When a task finishes, its bar stays where it is and the next task takes the
line below, so the block grows downward in the order the tasks ran, with the
step of the workflow on top.
A finished bar shows the time its task took, instead of the time still to
come, which also tells it apart from a task still running::

    cluster               |█████████     | step 5/8      [01:04<00:38]
    aligning reads        |██████████████| sample 2/2    [00:11]
    identifying mutations |██████████████| reference 2/2 [00:03]
    filtering regions     |██████████████| region 2/2    [00:11]
    clustering            |███████       | region 1/2    [00:18<00:18]

The block may not fill the terminal, or there would be no room for the log;
once it is full, each further finished bar goes into the log instead, so that
a long run does not push everything else off the screen.

Progress is tracked down to the level of one sample (or one file), not below
it.
The parts of a single sample's analysis, such as batches of reads or runs of
the clustering algorithm, begin and end too quickly to follow, so they get no
bar of their own.

Log messages appear above the bars, which stay at the bottom of the terminal
however many messages are written, even with ``-vv`` and even while tasks are
running in parallel.
When the command finishes, every bar stays on the terminal as a record of the
run: what each task did and how long it took.

Progress bars are hidden automatically when they would not be useful:

- with ``--quiet`` (``-q``) or quieter, whose purpose is to write little or
  nothing to the terminal
- when standard error is not a terminal, e.g. redirected to a file or piped
  into another program
- inside the Python API, which has no terminal of its own

Use ``--no-progress`` (between ``seismic`` and the subcommand) to hide them in
every case.


Log files
--------------------------------------------------------------------------------

By default, messages are also written to a log file at::

    ./log/seismic-rna_YYYY-MM-DD_hh-mm-ss.log

Use ``--log PATH`` (between ``seismic`` and the subcommand) to choose a
different path, or ``--log ""`` to disable log files.

The log file records messages at **the same verbosity level as the terminal**,
so ``-v`` enlarges both the terminal output and the log file.
Each line in the log file is prefixed with a timestamp, followed by the level,
the ID of the process that wrote the message, and the message itself::

    2026-06-04 14:22:01.123456  Info      12345 Clustering sars2-fse / region fse

The process ID is helpful for telling apart messages from parallel work: when
SEISMIC-RNA uses multiple CPUs at once, each one runs in its own process, so
messages that share a process ID come from the same task.

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
