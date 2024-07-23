
Clean FASTA Files
------------------------------------------------------------------------

The FASTA file cleaner modifies FASTA files so that they are compatible
with SEISMIC-RNA.

Background about cleaning FASTA files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Purpose of cleaning FASTA files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Variants of FASTA format exist that are not compatible with SEISMIC-RNA.
FASTA files from databases such as `NCBI`_ and `Gencode`_ often encode
metadata in their header lines, which SEISMIC-RNA cannot handle.
Additionally, some FASTA files may contain characters besides the four
nucleotides and ``N``, such as the `IUPAC extended nucleotide codes`_.
For details on the valid FASTA format for SEISMIC-RNA and the rationale
behind it, see :doc:`../../formats/data/fasta`.

Algorithm for cleaning FASTA files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The FASTA file cleaner of SEISMIC-RNA follows the following procedure:

- If a line begins with ``>``, then it is assumed to be a header line;
  otherwise, a body (i.e. sequence-encoding) line.
- If the line is a header line, then it is trimmed from the right side
  until all characters in the line (except the initial ``>`` character)
  are valid path characters.
- If the line is a body line, then all whitespace characters are removed
  and the remaining characters are converted to uppercase.
  If the FASTA sequence type is DNA, then every ``U`` is replaced with
  ``T``; and if RNA, then every ``T`` is replaced with ``U``.
  All remaining illegal characters that belong to the set of standard
  `IUPAC extended nucleotide codes`_ are replaced with ``N``.
  If any illegal characters remain after that, then an error is raised.

How to clean FASTA files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Command line for cleaning FASTA files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Clean the FASTA file ``refs.fa``::

    seismic +cleanfa refs.fa

This command will output a cleaned FASTA file with the same name as the
input FASTA file, located in the output directory (default ``./out``).
In this case: ``out/refs.fa``.

To change the output directory, use the ``--out-dir/-o`` option::

    seismic +cleanfa -o clean refs.fa

will output the cleaned FASTA file ``clean/refs.fa``.

By default, this command will not overwrite existing FASTA files.
To force overwriting, use the option ``--force``::

    seismic +cleanfa --force refs.fa

Troubleshooting cleaning FASTA files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Errors while cleaning FASTA files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

*Line contains sequence characters other than whitespace or IUPAC codes*

This error means that the sequence contains characters that cannot be
cleaned because they do not belong in FASTA files, so it is unclear how
they should be handled automatically.
The only characters that can be cleaned automatically are as follows:

- whitespace: deleted
- lowercase: converted to uppercase
- ``U`` (for DNA): replaced with ``T``
- ``T`` (for RNA): replaced with ``U``
- ``BDHKMRSVWY`` (`IUPAC extended nucleotide codes`_): replaced with ``N``
- standard four nucleotides, plus ``N``: unchanged

If the FASTA file contains any other characters, then first remove them
manually and then rerun ``+cleanfa``.
The error message will indicate which line(s) have invalid characters.


.. _`NCBI`: https://www.ncbi.nlm.nih.gov/
.. _`Gencode`: https://www.gencodegenes.org/
.. _`IUPAC extended nucleotide codes`: https://www.bioinformatics.org/sms/iupac.html
