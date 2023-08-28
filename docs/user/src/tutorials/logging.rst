
Messages and Log Files
========================================================================

As SEISMIC-RNA runs, it outputs messages about its normal activities and
any warnings or errors that occur. All messages are saved to a log file
to provide a record of what happened during the run. The most salient
messages are also printed to `standard output`_. This section explains
how to choose which messages to print and where to write log files.

Levels of messages
------------------------------------------------------------------------

Messages range from critical errors that halt the program to meticulous
details of primary interest to software developers. SEISMIC-RNA assigns
each message one of five levels from `Python's built-in logging tools`_:

1.  DEBUG: Detail of a routine event, e.g. the result of a calculation.
2.  INFO: Normal event or update on the status of the program.
3.  WARNING: Abnormal event or minor problem that causes no failure.
4.  ERROR: Problem causing one part of the program to fail.
5.  CRITICAL: Severe problem forcing the entire program to stop.

Controlling messages on standard output
------------------------------------------------------------------------



All levels of messages are written to the log file, which can be named using the ``--log`` (CLI) or ``log=`` (API) arguments (default: ``dreem_YYYY-MM-DD_hh-mm-ss.log``).
The level of messages printed to standard output can be controlled with the verbose and quiet arguments:

- Double Verbose (``-vv`` or ``--verbose --verbose``): all messages
- Verbose (``-v`` or ``--verbose``): all except debug messages
- Default: all except info and debug messages
- Quiet (``--quiet``): only errors and critical errors
- Double quiet (``--quiet --quiet``): only critical errors (not recommended)



.. _standard output: https://en.wikipedia.org/wiki/Standard_streams#Standard_output_(stdout)
.. _Python's built-in logging tools: https://docs.python.org/3/howto/logging.html
