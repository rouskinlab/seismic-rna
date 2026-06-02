********************************************************************************
seismic cleanfa
********************************************************************************


Purpose
================================================================================

``seismic cleanfa`` fixes FASTA files that are incompatible with SEISMIC-RNA.
Common issues from databases like NCBI and Gencode include metadata in header
lines and non-standard characters in sequences.
It trims headers to valid path characters, uppercases sequence bases, replaces
ambiguous IUPAC codes with ``N``, and converts ``U``↔``T`` as needed.


Inputs
================================================================================

FASTA files
    One or more FASTA files to clean.


Outputs
================================================================================

Cleaned FASTA files with the same names as the inputs, written to the output
directory (default ``./out``).


Quick example
================================================================================

Clean a reference FASTA::

    seismic cleanfa refs.fa

Output goes to ``out/refs.fa``.
Specify a different directory with ``-o``::

    seismic cleanfa -o clean refs.fa


Options
================================================================================

``--inplace/--no-inplace``
    Overwrite the original files instead of writing to ``--out-dir`` (default off).

    .. warning::

        ``--inplace`` is irreversible — back up your files first.

``--out-dir DIR`` (``-o``)
    Write cleaned files here (default ``./out``).
    Ignored when ``--inplace`` is set.

``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

*Line contains sequence characters other than whitespace or IUPAC codes*
    The sequence has characters that cannot be interpreted automatically.
    Remove them manually and rerun.
    The error message names the affected line(s).


See also
================================================================================

- :doc:`/formats/data/fasta` — the FASTA format required by SEISMIC-RNA
- :doc:`/use/workflow/align` — uses the cleaned FASTA as input


.. _NCBI: https://www.ncbi.nlm.nih.gov/
.. _Gencode: https://www.gencodegenes.org/
