
Logging
=======

Levels of Logging
-----------------

DREEM automatically logs messages of varying levels of importance to both standard output and to a log file.
The levels of importance match the levels defined by `Python's logging utilities <https://docs.python.org/3/howto/logging.html>`_:

1. CRITICAL: A file was missing, corrupted, misformatted, or otherwise unreadable.
2. ERROR: An item in a file (e.g. a read in a FASTQ file) could not be processed.
3. WARNING: A request from the user could not be fulfilled, but the final results will not be affected.
4. INFO: A step in the program (typically non-trivial or with file I/O) ran as expected.
5. DEBUG: A detail about a step or process (typically trivial and without file I/O).

Controlling logging output
--------------------------

All levels of messages are written to the log file, which can be named using the ``--log`` (CLI) or ``log=`` (API) arguments (default: ``dreem_YYYY-MM-DD_hh-mm-ss.log``).
The level of messages printed to standard output can be controlled with the verbose and quiet arguments:

- Double Verbose (``-vv`` or ``--verbose --verbose``): all messages
- Verbose (``-v`` or ``--verbose``): all except debug messages
- Default: all except info and debug messages
- Quiet (``--quiet``): only errors and critical errors
- Double quiet (``--quiet --quiet``): only critical errors (not recommended)
