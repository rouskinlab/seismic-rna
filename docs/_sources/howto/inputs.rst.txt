
List Input Files
========================================================================

For commands that take a list of input files as positional arguments,
these files can be given in three ways -- or combinations thereof --
so that you can choose a convenient method to list the input files.

How to list input files
------------------------------------------------------------------------

- If you have a small number of input files, then the easiest method is
  to list every input file explicitly::

    seismic {command} {file-1} {file-2} {file-3}

- If you have a large number of input files that have similar names,
  then you can use `glob patterns`_ to list all of them.
  For example, if you have 83 files named ``file-1``, ``file-2``, and
  so on up to ``file-83``, then you can process all of them with ::

    seismic {command} file-*

  Note that glob pattern matching is a general ability of the shell, not
  a special feature of SEISMIC-RNA.
  The shell itself searches for all files that match the glob pattern,
  then implicitly replaces the pattern with a space-separated list of
  files that match; SEISMIC-RNA never "sees" the original glob pattern.
  Refer to `glob patterns`_ for more information.

- If your input files are within a directory, then you can give the path
  of the directory, which will be searched for files recursively, with
  no limit to the depth of the search.
  This method is particularly useful if you want to process all files
  in a directory, such as your output directory (assume it is ``out``)::

    seismic {command} out

You can also combine any of the above methods at your convenience.
For example, to process all files in the directories ``out/sample-1``,
``out/sample-2``, and ``out/sample-3``, as well as ``out/sample-26``;
plus the files ``out/sample-7/mask/ref-6/full/mask-report.json``,
``out/sample-7/cluster/ref-6/full/cluster-report.json`` and
``out/sample-9/relate/ref-3/relate-report.json``, you could use::

    seismic {command} out/sample-[123] out/sample-26 out/sample-7/*/ref-6/full/*-report.json out/sample-9/relate/ref-3/relate-report.json


.. note::
    `Glob patterns`_ work *only* with positional arguments, not optional
    arguments (i.e. those given via option such as ``--fastqx``).
    Thus, you *cannot* use glob patterns to list FASTQ files via the
    options ``-x``, ``-y``, ``-z``, ``-X``, ``-Y``, and ``-Z``.
    You can, however, pass both individual FASTQ files and directories
    (to be searched recursively for FASTQ files) to these options.

.. _glob patterns: https://en.wikipedia.org/wiki/Glob_(programming)
