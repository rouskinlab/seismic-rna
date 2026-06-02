
Specify Input Files
================================================================================

Commands that accept input files as positional arguments support three ways to
list them, which you can combine freely.


List files explicitly
--------------------------------------------------------------------------------

The simplest method: name each file directly::

    seismic {command} out/sample-1/idmut/ref-1 out/sample-2/idmut/ref-1


Use glob patterns
--------------------------------------------------------------------------------

Use shell wildcards to match many files at once::

    seismic {command} out/*/idmut/ref-1

The shell expands the pattern before passing the list to SEISMIC-RNA.
See `glob patterns`_ for the full syntax (``*`` matches any text, ``?``
matches one character, ``[12]`` matches either ``1`` or ``2``).


Pass a directory
--------------------------------------------------------------------------------

Passing a directory causes SEISMIC-RNA to search it recursively for all
matching files::

    seismic {command} out/

This is useful when you want to process every output file from a previous step.


Combine methods
--------------------------------------------------------------------------------

You can mix all three methods in one command::

    seismic {command} out/sample-[123] out/sample-26 out/sample-7/filter/ref-1


Glob patterns in options that accept multiple paths
--------------------------------------------------------------------------------

Several options that accept more than one file or directory can also
expand glob patterns. The options that support this are:

- ``-x`` / ``--fastqx``, ``-y`` / ``--fastqy``, ``-z`` / ``--fastqz``
  (FASTQ inputs)
- ``-X`` / ``--dmfastqx``, ``-Y`` / ``--dmfastqy``,
  ``-Z`` / ``--dmfastqz`` (demultiplexed FASTQ inputs)
- ``-d`` / ``--param-dir`` (simulation parameter directories)
- ``--ct-file``, ``--ct-pos-5`` (connectivity table files)
- ``--struct-file`` (reference structure files)
- ``--mask-pos-file``, ``--drop-read-file`` (filter input lists)

.. important::

    You **must prefix every wildcard character with a backslash** (``\``)
    when the pattern is the value of an option flag. Without the
    backslash, the shell expands the pattern into separate arguments
    *before* SEISMIC-RNA sees them: only the first matched file becomes
    the value of the option, and the rest are left over as stray
    arguments, which causes the command to fail.

For example, to pass both mates of a paired-end sample as
``-x``, escape the brackets so SEISMIC-RNA (not the shell) does the
expansion::

    seismic align refs.fa -x samples/sample1_R\[12\].fastq

The same rule applies to ``*`` and ``?``::

    seismic align refs.fa -x samples/sample1_R\*.fastq

You can repeat the option to combine literal paths and escaped
patterns::

    seismic align refs.fa -x sample1_R1.fastq -x sample1_R2.fastq \
                         -x extras/\*.fastq

A pattern that matches no files raises an error, so typos in patterns
are surfaced immediately rather than silently ignored.

Options that accept only a single file (for example ``--regions-file``,
``--samples-meta``, ``--refs-meta``) do **not** expand glob
patterns; give them an exact path.


.. _glob patterns: https://en.wikipedia.org/wiki/Glob_(programming)
