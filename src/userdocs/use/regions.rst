Define Regions
================================================================================

A region is a contiguous range of positions within one reference sequence.
Most analysis steps can restrict their work to one or more such regions rather
than the full reference length.


Define a region by coordinates
--------------------------------------------------------------------------------

Specify the first (5') and last (3') position using 1-based numbering (the
first base of the reference is position 1).
Both endpoints are included.

Use ``--region-coords`` (``-c``) followed by the reference name and the two
coordinates::

    seismic filter -c ref-1 34 71 out/

defines a region spanning positions 34–71 of ``ref-1``.
Repeat ``-c`` to define multiple regions.


Define a region by primer sequences
--------------------------------------------------------------------------------

For amplicon samples prepared with RT-PCR, it is usually more convenient to
define a region by primer sequences than by coordinates.
SEISMIC-RNA locates the exact primer-binding sites in the reference (no
mismatches are permitted) and creates a region spanning from one base
downstream of the forward primer's 3' end to one base upstream of the reverse
primer's 5' end.

Use ``--region-primers`` (``-P``) followed by the reference name and the two
primer sequences (give the reverse primer as written on the oligonucleotide
itself, not as its reverse complement)::

    seismic filter -P ref-1 TCGAAT GTTACG out/

To exclude a few extra bases near each primer end (to reduce edge artifacts),
add ``--primer-gap N`` (default 0).


Define regions in a file
--------------------------------------------------------------------------------

To define many regions reproducibly or share definitions across multiple runs,
list them in a CSV file and pass it with ``--regions-file`` (``-i``)::

    seismic filter -i regions.csv out/

See :doc:`/formats/meta/regions` for the CSV format.


Multiple regions per reference
--------------------------------------------------------------------------------

Repeat ``-c``, ``-P``, or rows in a ``--regions-file`` to define as many
regions as you need for each reference.
Regions may overlap.
Regions for one reference have no effect on any other reference.


Region names
--------------------------------------------------------------------------------

- Regions from ``-c`` or ``-P`` are named ``{first}-{last}`` (e.g. ``34-71``).
- Regions from ``--regions-file`` use the name in the ``Region`` column.
- If no region is defined for a reference, one region named ``full`` spanning
  the entire reference is created automatically.


See also
--------------------------------------------------------------------------------

- :doc:`/formats/meta/regions` — CSV format for ``--regions-file``
- :doc:`inputs` — ways to pass multiple input files to a command
