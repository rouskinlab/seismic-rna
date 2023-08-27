
Relation Vectors
------------------------------------------------------------------------

Primary relationships of reads and references
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A relation vector encodes the relationship between one sequencing read
(or pair of mated reads) and each base in the reference sequence.
SEISMIC-RNA defines eight primary relationships:

- **Match**: The reference base aligned to a high-quality base of the same kind in the read.
- **Deletion**: The reference base aligned to the space between two bases in the read.
- **5' of Insertion**: The reference base aligned to any read base immediately 5' of an extra base in the read.
- **3' of Insertion**: The reference base aligned to any read base immediately 3' of an extra base in the read.
- **Substitution to A**: The reference base is not A and aligned to a high-quality A in the read.
- **Substitution to C**: The reference base is not C and aligned to a high-quality C in the read.
- **Substitution to G**: The reference base is not G and aligned to a high-quality G in the read.
- **Substitution to T**: The reference base is not T and aligned to a high-quality T in the read.

This figure illustrates six of these primary relationships, as well the
"blank" relationship for positions in the reference that lie outside the
span of the read.

.. image::
    relationships.png

.. _encode-primary-relation

Encoding primary relationships
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each position in a relation vector is physically represented by one byte
(eight bits); each bit corresponds to one of the eight types of primary
relationship. For a given position in the relation vector, the byte(s)
corresponding to the relationship(s) at that position are turned on (set
to 1); all other bits are set to 0. The following table indicates which
relationship each bit represents. Each bit is shown as the sole 1 within
an entire byte, an eight-digit binary number. The number's decimal (Dec)
and hexadecimal (Hex) forms are also shown.

========== ===== ===== ===================
 Byte       Dec   Hex   Relationship
========== ===== ===== ===================
 00000001   001    01   Match
 00000010   002    02   Deletion
 00000100   004    04   5' of Insertion
 00001000   008    08   3' of Insertion
 00010000   016    10   Substitution to A
 00100000   032    20   Substitution to C
 01000000   064    40   Substitution to G
 10000000   128    80   Substitution to T
========== ===== ===== ===================

Ambiguous relationships
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Oh, if only encoding relation vectors were that straightforward! Though
most bases in a read will have primary relationships with the bases in
the reference to which they align, two phenomena make it more difficult
for some bases to define the relationship:

- low-quality base calls
- ambiguous insertions and deletions

Low-quality base calls
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A low-quality base call is defined as having a `Phred quality score`_
below the user-specified threshold (default: 25). Low-quality base calls
are processed as if they could be any of the four nucleotides. The type
of relationship that would occur if the read base were each of the four
nucleotides is determined, and then the results are summed. For example,
suppose a low-quality base call aligns to a T in the reference. If the
read base (which we don't know, because it is low-quality) were actually
an A, then the relationship would be a substitution to A. Likewise for
C and G. If the read base were actually a T, then the relationship would
be a match. So the four possible primary relationships are summed:

- Sub. A (``00010000``)
- Sub. C (``00100000``)
- Sub. G (``01000000``)
- Match (``00000001``)

The sum is ``01110001``; this becomes the byte in the relation vector.

.. note::
    More accurately, the operation is `bitwise or`_, not addition. These
    two operations produce identical results for primary relationships
    because none of the eight bits are set to 1 in more than one primary
    relationship byte (see :ref:`encode-primary-relation`_).

Each row of the following table repeats this calculation for a different
reference base (Ref). Each column "Read: A/C/G/T?" shows what would be
the relationship if the low-quality base in the read were actually the
base in the column header. For example, in the first row, the reference
base is A: if the read base were A (first column), then the relationship
would be a match (``00000001``); and if it were C (second column), then
the relationship would be a substitution to C (``00100000``). The column
"Byte" shows the resulting ambiguous relationship, obtained by summing
the four columns "Read: A/C/G/T?". Its decimal (Dec) and hexadecimal
(Hex) values are also shown.

===== ========== ========== ========== ========== ========== ===== =====
 Ref   Read: A?   Read: C?   Read: G?   Read: T?   Byte       Dec   Hex
===== ========== ========== ========== ========== ========== ===== =====
  A    00000001   00100000   01000000   10000000   11100001   225    e1
  C    00010000   00000001   01000000   10000000   11010001   209    d1
  G    00010000   00100000   00000001   10000000   10110001   177    b1
  T    00010000   00100000   01000000   00000001   01110001   113    71
===== ========== ========== ========== ========== ========== ===== =====

.. note::
    A byte that has more than one bit set to 1 does **not** count more
    than once towards the total number of matches or mutations. To learn
    how mutations in relation vectors are counted, see [REF].

Ambiguous insertions and deletions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Insertions and deletions (collectively, "indels") in the read can cause
ambiguities that even the highest quality sequencing reads could not
prevent. When one or more bases are inserted or deleted in a repetitive
sequence, the exact base that mutated cannot be determined. For example,
if the reference is ``ACTTGA`` and the read is ``ACTGA``, then one T was
clearly deleted from the read, but determining whether it was the first
or second T is impossible. The alignments are equally good:

**Deletion of the first T**::

    AC-TGA
    || |||
    ACTTGA

**Deletion of the second T**::

    ACT-GA
    ||| ||
    ACTTGA

Ambiguities in the location of a relationship are encoded by turning on
the bit of every possible relationship at each position. In the above
example, position 3 of the reference could be a deletion (``00000010``)
or a match (``00000001``), so the byte it receives is the `bitwise or`_
of the two relationships: ``00000011``. Likewise for position 4. Thus,
the relationship vector for the above alignment would be

===== ========== =====
 Pos   Byte       Hex
===== ========== =====
  1    00000001    01
  2    00000001    01
  3    00000011    03
  4    00000011    03
  5    00000001    01
  6    00000001    01
===== ========== =====


.. _Phred quality score: https://en.wikipedia.org/wiki/Phred_quality_score
.. _bitwise or: https://en.wikipedia.org/wiki/Bitwise_operation#OR
