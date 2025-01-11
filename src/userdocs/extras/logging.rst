
Log Messages
================================================================================

As SEISMIC-RNA runs, it logs messages about routine events and problems.
All messages are saved to a log file to provide a record of the run.
You can control which of these messages are also printed to `standard error`_,
as described in this section.

Logged messages come in eight levels
--------------------------------------------------------------------------------

Messages range from fatal errors to details of primary use in troubleshooting:

======= ========== =============================================================
 Level   Name       Used for
======= ========== =============================================================
 4       DETAIL     Detail about a normal event
 3       ROUTINE    Normal internal event (e.g. regular function call)
 2       ACTION     Normal external event (e.g. file output or shell command)
 1       TASK       Task that can be run in series or parallel
 0       STATUS     Current command or step of the workflow
 -1      WARNING    Abnormal event or problem that can be fully recovered
 -2      ERROR      Problem preventing a requested output from being written
 -3      FATAL      Problem preventing all requested outputs from being written
======= ========== =============================================================

Logging messages on standard error
--------------------------------------------------------------------------------

Format of logged messages on standard error
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each message on standard error consists of the message level, in all capital
letters, a tab character, and the message itself, for example::

    DETAIL  This text shows the format of a message on standard error.

Levels are also logged in different colors to make them easy to distinguish.

How to control which messages are logged to standard error
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the only messages that are logged to standard error are of level 0
or more negative (fatal, error, warning, and command).
The flags ``--verbose`` (``-v``) and ``--quiet`` (``-q``) change the threshold
level for logging to standard error:

======= =========== ============================================================== =============================
 Flag    Threshold   Logs to standard error                                          Useful for
======= =========== ============================================================== =============================
 -vvvv   4           fatal, error, warning, status, task, action, routine, detail   Detailed troubleshooting
 -vvv    3           fatal, error, warning, status, task, action, routine           General troubleshooting
 -vv     2           fatal, error, warning, status, task, action                    Tracking file operations
 -v      1           fatal, error, warning, status, task                            Monitoring individual tasks
         0           fatal, error, warning, status                                  Monitoring overall status
 -q      -1          fatal, error, warning                                          Focusing on all problems
 -qq     -2          fatal, error                                                   Focusing on errors
 -qqq    -3          fatal                                                          Focusing on fatal errors
 -qqqq   -4                                                                         Silencing all logging
======= =========== ============================================================== =============================

The flags must come between the command ``seismic`` and the sub-command.
For example, to run the ``cluster`` step with double-verbose logging::

    seismic -vv cluster -k 3 out/sars2-fse

This is because ``--verbose`` and ``--quiet`` affect the global logging
system; they are not interpreted separately by each sub-command.
The one exception is the ``test`` command, which uses the `Python unittest`_
framework that has its own, separate system for controlling verbosity, e.g. ::

    seismic test -vv

Logging messages to files
--------------------------------------------------------------------------------

All messages are written to the log file.

Format of messages in log files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each message in a log file consists of the label ``LOGMSG>``, the message level,
the date (YYYY-MM-DD) and the time (hh:mm:ss.microseconds).
On the next line(s) is the message itself, after which comes one blank line::

    LOGMSG> WARNING on 2024-04-08 at 15:20:09.277928
    It appears that something has blotted out the sun in Rochester, NY.

How to control the path to log files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default path of the log file is ::

    ./log/seismic-rna_YYYY-MM-DD_hh-mm-ss.log

where the date and time are the time at which the command began running.
A different path can be given with the ``--log`` option.
Leaving it blank (``--log ""`` or ``--log ''``) disables writing log files.
As with the ``--verbose`` and ``--quiet`` options, the ``--log`` option is
global and must be given between ``seismic`` and the sub-command::

    seismic --log log/final-2_use-this-one.log cluster -k 3 out/sars2-fse

How to view log files in real time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regardless of the verbosity or quietness specified by the ``-v``/``-q`` flags,
all messages from all levels can be viewed in real time within the log file.
To do so, open the log file in the text viewer `less`_::

    less log/seismic-rna_2024-04-08_15-20-09.log

Press Shift-F; this message should appear at the bottom of the screen::

    Waiting for data... (interrupt to abort)

As new messages are written to the file, they will appear immediately at the
bottom of the screen.
To stop updating the file, press Ctrl-C.
Then, to quit out of less, press the Q key.


.. _standard error: https://en.wikipedia.org/wiki/Standard_streams#Standard_error_(stderr)
.. _Python's built-in logging tools: https://docs.python.org/3/howto/logging.html
.. _Python unittest: https://docs.python.org/3/library/unittest.html
.. _less: https://greenwoodsoftware.com/less/
