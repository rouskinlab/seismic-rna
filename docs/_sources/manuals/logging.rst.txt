
Log Messages
========================================================================

As SEISMIC-RNA runs, it outputs messages about its normal activities and
any warnings or errors that occur. All messages are saved to a log file
to provide a record of what happened during the run. The most salient
messages are also printed to `standard output`_. This section shows how
to control which messages are printed and where log files are written.

Logged messages come in five levels
------------------------------------------------------------------------

Messages range from critical errors that halt the program to meticulous
details of primary interest to software developers. SEISMIC-RNA assigns
each message one of five levels from `Python's built-in logging tools`_:

======= ========== =====================================================
 Level   Name       Used for
======= ========== =====================================================
 1       DEBUG      Detail about a routine event
 2       INFO       Normal event or update on the status of the program
 3       WARNING    Abnormal event or minor problem causing no failure
 4       ERROR      Problem causing one part of the program to fail
 5       CRITICAL   Severe problem forcing the entire program to stop
======= ========== =====================================================

Logging messages on standard output
------------------------------------------------------------------------

Format of logged messages on standard output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each message on standard output consists of the message level, in all
capital letters, a tab character, and the message itself, for example::

    INFO    This text shows the format of a message on standard output.

Each level is also output in a different color to make it easier to spot
more important messages.

How to control which messages are logged to standard output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the only messages that are printed to standard output are of
level 3 or higher (warning, error, and critical). This means that, if no
problems occur, then SEISMIC-RNA will print nothing. The ``--verbose``
(``-v``) and ``--quiet`` (``-q``) flags, respectively, lower and raise
the threshold level for printing messages:

====== =========== ======================================= =================================
 Flag   Threshold   Prints                                  Useful for
====== =========== ======================================= =================================
 -vv    1           critical, error, warning, info, debug   Troubleshooting  details
 -v     2           critical, error, warning, info          Monitoring status and progress
 none   3           critical, error, warning                Noticing all problems (default)
 -q     4           critical, error                         Ignoring minor problems
 -qq    5           critical                                Ignoring all non-fatal errors
====== =========== ======================================= =================================

The flags must come between the command ``seismic`` and the sub-command.
For example, to run the ``cluster`` step with double-verbose logging::

    seismic -vv cluster -k 3 out/sars2-fse

This is because ``--verbose`` and ``--quiet`` affect the global logging
system; they are not interpreted separately by each sub-command. The one
exception is the ``test`` sub-command, which uses the `Python unittest`_
framework that has its own, separate system for controlling verbosity::

    seismic test -vv

Logging messages to files
------------------------------------------------------------------------

All messages are written to the log file.

Format of messages in log files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each message in a log file consists of the label ``LOGMSG>``, the date
and time (YYYY-MM-DD hh:mm:ss,dcm), the module that printed the message,
and the message level, separated by tab characters on one line. On the
next line(s) is the message itself, after which comes one blank line::

    LOGMSG> 2024-04-08 15:20:09,277 seismicrna.core.shell   WARNING
    It appears that something has blotted out the sun in Rochester, NY.

How to control the path to log files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default path of the log file is ::

    ./log/seismic-rna_YYYY-MM-DD_hh-mm-ss.log

where the date and time are the time at which the command began running.
A different path can be given with the ``--log`` option; ``--log ""`` or
``--log ''`` disables writing log files. As with the ``--verbose`` and
``--quiet`` options, the ``--log`` option is global and must be given
between ``seismic`` and the sub-command::

    seismic --log log/final-2_use-this-one.log cluster -k 3 out/sars2-fse

How to view log files in real time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Regardless of the verbosity or quietness specified by the ``-v``/``-q``
flags, all messages from all levels can be viewed in real time within
the log file. To do so, open the log file in the text viewer `less`_::

    less log/seismic-rna_2024-04-08_15-20-09.log

Press Shift-F; this message should appear at the bottom of the screen::

    Waiting for data... (interrupt to abort)

As new messages are written to the file, they will appear immediately
at the bottom of the screen. To stop updating the file, press Ctrl-C.
Then, to quit out of less, press the Q key.


.. _standard output: https://en.wikipedia.org/wiki/Standard_streams#Standard_output_(stdout)
.. _Python's built-in logging tools: https://docs.python.org/3/howto/logging.html
.. _Python unittest: https://docs.python.org/3/library/unittest.html
.. _less: https://greenwoodsoftware.com/less/
