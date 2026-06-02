
Algorithm for Finding Ambiguous Insertions and Deletions
========================================================================

Background on Ambiguous Insertions and Deletions
------------------------------------------------------------------------

When a read is aligned to a reference sequence, the aligner sometimes
must insert or delete one or more bases -- collectively, "indels" -- to
make the read line up with the reference.
In a region where the same base (or pattern of bases) repeats, the place
where the aligner puts an indel is often arbitrary, because several
different placements describe the read equally well.

For example, consider a reference that contains a run of four ``A`` bases
and a read in which one of those ``A`` bases is missing::

    Position   1 2 3 4 5 6
    Reference  C A A A A G
    Read       C A A A   G

The single missing ``A`` could be any of the four ``A`` bases in the
reference: deleting the first, second, third, or fourth ``A`` produces
exactly the same read.
The deletion is therefore *ambiguous* -- it has four equally valid
positions.
The aligner reports only one of them (chosen arbitrarily), even though
all four fit the data equally well.

The same situation arises for insertions.
If the read instead contained five ``A`` bases where the reference has
four, the extra ``A`` could be placed at any point within the run.

Why this matters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SEISMIC-RNA counts mutations -- including indels -- at each position in
order to measure the chemical probing signal along the molecule.
If an ambiguous indel were assigned to only the single position that the
aligner happened to choose, then that one position would appear to carry
a definite mutation, while the neighboring positions that are equally
likely would appear unmutated.
This false precision would distort the signal.

To avoid it, SEISMIC-RNA finds every position at which each indel could
equivalently be placed and marks all of them as ambiguous.
Downstream steps can then treat these positions as uncertain rather than
trusting the aligner's arbitrary choice.

Goal of the algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given the indels that the aligner placed in a read, find every reference
position at which each indel could be placed *without changing how well
the read matches the reference*, and mark all of those positions as
possible indels.

How the Algorithm Works
------------------------------------------------------------------------

Sliding an indel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The core operation is to *slide* an indel by one position, swapping it
with the neighboring base.
A deletion slides along the reference (the deleted base is a reference
base); an insertion slides along the read (the inserted base is a read
base).

An indel may slide by one position only if the resulting alignment fits
the read equally well.
In practice, this means that the base the indel trades places with must
relate to the reference the same way both before and after the move: a
match must stay a match, and a substitution must stay a substitution.
If sliding would turn a match into a mismatch -- or otherwise change how
a base relates to the reference -- then the two placements are *not*
equivalent, so the indel cannot slide past that point.
In this way, an indel moves freely within a repeat but stops at its
boundaries.

Bases that were read at low quality, or that are ambiguous (``N``), are
treated as compatible with any reference base, so an indel can slide
across them.
An uncertain base call should not artificially limit the range of
ambiguity.

Keeping indels of different kinds separate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before any sliding, the indels in a read are organized into groups, where
each group holds only indels of a single kind -- all insertions or all
deletions -- with no indel of the other kind between them.
This grouping serves two purposes.
First, it keeps insertions and deletions from merging into one another as
they slide, since an insertion and a deletion are distinct events that
must remain separate.
Second, when a repeat contains more than one indel of the same kind, the
members of a group are allowed to slide *through* one another; handling
them together keeps track of their relative order as they pass.

The two sliding passes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The algorithm maps out the full range of each indel's ambiguity in two
passes that move in opposite directions.

The **first pass** slides every indel as far toward the 5' end as it can
go.
The purpose of this pass is simply to bring each indel to one extreme --
the 5'-most position it could possibly occupy -- so that the second pass
has a complete, well-defined starting point.
This pass works through the groups from the most 5' to the most 3', so
that a group still moving 5' is never blocked by a group behind it that
has not yet moved out of the way.

The **second pass** then slides every indel step by step toward the 3'
end, as far as it can go, and records each position it occupies along the
way.
Because each indel begins this pass at its 5'-most position and travels to
its 3'-most position, the second pass visits -- and marks -- every
position at which the indel could equivalently be placed.
This pass works through the groups in the opposite order, from the most
3' to the most 5', so that an indel moving 3' is not blocked by another
indel that has not yet moved out of its path.

Checks made before each move
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each time the algorithm tries to slide an indel by one position, it makes
three checks, and performs the move only if all three pass:

- **Boundary:** The indel must not move outside the part of the reference
  (for a deletion) or the read (for an insertion) that the read actually
  covers.
  An indel cannot be placed beyond the ends of the alignment.
- **Collision:** The indel must not move into or next to a position held
  by an indel of the opposite kind in a neighboring group.
  A deletion and an insertion must never overlap, as that would merge two
  distinct events into one.
- **Consistency:** The move must leave the read fitting the reference
  equally well, as described above.
  This is the check that actually defines the limits of the ambiguity.

When an indel does move, both positions involved are annotated: the
position it leaves is marked with what the read implies there if no indel
were present (for instance, a match or a substitution), and the position
it enters is marked as a possible indel.
Because the second pass carries each indel across its whole range, every
position the indel could occupy ends up marked as a possible indel.

A note on insertions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An inserted base lies *between* two reference positions rather than at one
of them.
SEISMIC-RNA therefore records each insertion at one of the two reference
positions that flank it; a setting determines whether the position on the
3' or the 5' side is used, and that choice is applied consistently to
every insertion.
Apart from this bookkeeping, an insertion slides within a repeat and is
bounded by the same kinds of checks as a deletion.

Example
------------------------------------------------------------------------

Suppose the reference contains a run of four ``A`` bases and a read is
missing one of them.
The aligner reports the deletion at the first ``A`` (reference position
2)::

    Position   1 2 3 4 5 6
    Reference  C A A A A G
    Read       C A A A   G
                 ^
                 deletion reported here

**First pass (slide 5').**
The deletion is already at position 2, the start of the run, so it cannot
move any further 5'.
Position 1 is ``C``, which the read matches, so sliding the deletion there
would turn a match into a mismatch -- an inconsistent move.

**Second pass (slide 3', recording each position).**
The algorithm now slides the deletion 3' one step at a time:

- Position 2 → 3: trades the deletion with an ``A``; an ``A`` matches an
  ``A`` either way, so the move is consistent.
- Position 3 → 4: consistent for the same reason.
- Position 4 → 5: consistent for the same reason.
- Position 5 → 6: the base just past the run is ``G``, which the read
  also reads as ``G`` (a match); moving the deletion onto it would break
  that match, so the move is inconsistent and the slide stops.

Every position the deletion occupied -- 2, 3, 4, and 5 -- is marked as a
possible deletion.
Each of these positions is also marked as a possible match (what it would
be if the deletion were elsewhere), so each carries the combined meaning
"either a deletion or a match here."
That combined annotation is exactly what it means for the position to be
ambiguous.

The result correctly captures that the missing base could be any one of
the four ``A`` bases in the run, rather than the single position the
aligner happened to report.
