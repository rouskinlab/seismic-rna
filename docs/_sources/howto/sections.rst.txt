
Define Sections
========================================================================

Reference sequences can be divided into (possibly overlapping) sections,
each of which is a contiguous range of positions in one sequence.

.. _sections_coords:

How to define sections using coordinates or primers
------------------------------------------------------------------------

How to define sections using coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A section can be defined by specifying its first and last coordinates,
with respect to the numbering of positions in the reference sequence.
SEISMIC-RNA uses the following conventions:

- The first (i.e. 5') position in the reference sequence is numbered 1,
  as is conventional for biological sequences.
  Consequently, the last (i.e. 3') position in the reference sequence
  is assigned the number equal to the length of the reference sequence.
- A section includes both the first and last coordinates.
  Note that the length of the section therefore equals the last minus
  the first coordinate plus one.
  For example, if the reference sequence is 10 nt long, then the section
  whose first/last coordinates are 4 and 8, respectively, will include
  positions 4, 5, 6, 7, and 8; but not positions 1, 2, 3, 9, or 10.
- The first coordinate must be a positive integer, and the last must be
  an integer greater than or equal to the first coordinate minus one and
  less than or equal to the length of the reference sequence.
- Setting the last coordinate to the first coordinate minus one creates
  a zero-length section; this behavior is permitted in order to handle
  this edge case smoothly, but rarely if ever has practical use.

On the command line, define a section using coordinates via the option
``--coords`` (``-c``) followed by the name of the reference and the two
coordinates.
For example, ``-c refA 34 71`` would define a section spanning positions
34 to 71 of reference "refA".

How to define sections using primer sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For samples that were prepared as amplicons (using `RT-PCR`_ with a pair
of primers flanking a specific sequence), it is usually easier to define
sections using the sequences of the forward and reverse primers compared
to using coordinates.
SEISMIC-RNA will compute the coordinates from the primer sequences using
the following conventions:

- The entire sequence of the forward primer and the reverse complement
  of the reverse primer must both occur exactly once in the reference
  sequence (no mismatches or gaps are permitted).
- Primers will cover up any mutations during `RT-PCR`_, so the sequences
  they bind provide no information for mutational profiling.
  Thus, SEISMIC-RNA makes the section start one position downstream of
  the 3' end of the forward primer, and end one position upstream of the
  5' end of the reverse complement of the reverse primer.
- Artifacts may occur near the ends of the primers.
  If so, then a number of positions at the both ends of the section can
  be ignored using the option ``--primer-gap {n}``, where ``{n}`` is the
  number of positions to ignore.

On the command line, define a section using primers via the option
``--primers`` (``-p``) followed by the name of the reference and the two
primer sequences; the reverse primer must be written as the sequence of
the oligonucleotide itself, not its reverse complement.
For example, if the sequence of "refA" is ``TTTCGCTATGTGTTAC``, then
``-p refA TCG AAC`` would define a section, depending on the primer gap,
spanning these positions:

- 6-11 (``CTATGT``) with ``--primer-gap 0`` (the default)
- 7-10 (``TATG``) with ``--primer-gap 1``
- 8-9 (``AT``) with ``--primer-gap 2``

How to define sections using a file of coordinates and/or primers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To create many sections at once, and to do so reproducibly, it is more
convenient to define them in a file than on the command line.
For information on this file, see :doc:`../formats/meta/sections`.
Provide a file defining sections using the option ``--sections-file``
(``-s``), for example ``-s sections.csv``.

How to define multiple sections simultaneously
------------------------------------------------------------------------

How to define multiple sections for one reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On the command line, the ``-c`` and ``-p`` options can be given multiple
times to define more than one section for a reference.
In a sections file (``-s``), each reference can be given multiple times,
as long as the name of each section of that reference is unique.

For example, typing ``-c refA 34 71 -c refA 56 89 -c refA 103 148`` on
the command would make three sections for reference "refA":

- positions 34 to 71
- positions 56 to 89
- positions 103 to 148

Note that sections are allowed to overlap or contain each other, as with
34-71 and 56-89.

How to define sections for multiple references
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sections are defined for each reference sequence separately and have no
effect on sections of any other reference sequence.
For example, suppose that for three references, "refA", "refB", and
"refC", one section is defined for "refA" and one for "refB" using the
option ``-c refA 34 71 -c refB 86 130``.

- "refA": positions 34 to 71
- "refB": positions 86 to 130

References "refB" and "refC" do not get section 34-71, nor do "refA" and
"refC" get section 86-130.
If no sections are defined for a reference sequence (here, "refC"), then
one section, spanning the full reference sequence and named "full", is
created automatically.
This applies to sections given as coordinates or primers on the command
line (with ``-c`` and ``-p``) as well as in a file (with ``-s``).

How to use the full length of a reference sequence automatically
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the full length of a reference sequence, simply avoid specifying
any coordinates or primers for that reference on the command line or in
a sections file.
Keep in mind that you may specify coordinates or primers for any other
references, since all references are sectioned separately.

How to name sections
------------------------------------------------------------------------

Sections are named as follows:

- If the section is defined in a file, then its name is taken from the
  "Section" column of the file (see :doc:`../formats/meta/sections`).
- If the section is defined on the command line by its coordinates or
  primers, then its name is the first and last coordinates, hyphenated
  (e.g. ``-c refA 34 71`` would create a section named "34-71").
- If the section is created automatically because no other sections were
  defined for its reference sequence, then its name is "full".

.. _RT-PCR: https://en.wikipedia.org/wiki/Reverse_transcription_polymerase_chain_reaction
