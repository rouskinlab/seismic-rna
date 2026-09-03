********************************************************************************
seismic duplex
********************************************************************************


Purpose
================================================================================

``seismic duplex`` combines two references into one **duplex** so that
:doc:`fold` can predict the structure the two molecules form *together*
(cofolding), not just the structure each one forms on its own.

Use it when the biological question involves two RNA molecules interacting,
for example:

- an mRNA and a small regulatory RNA that base-pair with each other
- a target RNA and an antisense oligonucleotide
- a homodimer, i.e. two copies of the same RNA pairing with each other

``seismic duplex`` itself does no structure prediction: it fuses the two
strands and their chemical probing data into a single profile, which
``seismic fold`` then cofolds with RNAcofold_.


How a duplex is built
--------------------------------------------------------------------------------

The two strands are joined into one sequence, 5' strand first, with a short
linker of three ``N`` bases between them::

    5' strand (refA)          linker   3' strand (refB)
    |-----------------------|  NNN   |-----------------------|
     positions 1-60           61-63   64-123

The linker only makes the duplex a single, continuously numbered molecule, so
that every downstream step (:doc:`fold`, :doc:`graph`, :doc:`draw`) can treat
it exactly like an ordinary folded region.
It carries no data, and it is removed before RNAcofold runs, so it never
affects the predicted structure; base pairs between the two strands simply
appear as a loop that spans the linker.

Each strand keeps its own mutational data, which end up on that strand's half
of the fused molecule.
The new duplex is named after both of its parts: reference ``refA__refB``
from references ``refA`` and ``refB``, and sample ``sample1__sample2`` if the
two strands came from different samples.
A duplex always spans its whole fused sequence, so its region is named
``full``.

If either strand is clustered, the duplex contains one profile for every
combination of the two strands' clusters (each strand at its best number of
clusters).
For example, a 2-cluster strand duplexed with a 3-cluster strand gives 6
duplex profiles, which are cofolded independently.


Inputs
================================================================================

Table CSV files or output directories
    Per-position mutation rate tables, the same inputs :doc:`fold` accepts
    (``filter-position-table.csv`` from :doc:`filter`, or
    ``cluster-position-table.csv`` from :doc:`cluster`).
    Pass a table CSV, the directory containing it, or a higher-level
    directory.  See :doc:`/use/inputs`.

    Duplex tables already present among the inputs are ignored, so re-running
    ``seismic duplex`` over a directory that holds its own output never fuses
    duplexes into larger chimeras.

A partner sequence without data (optional)
    The 3' strand does not need probing data of its own.
    It can instead be given as a FASTA file (``--duplex-file``) or typed
    directly on the command line (``--duplex-sequence``), which is the usual
    way to model an antisense oligonucleotide or any unprobed partner.


Outputs
================================================================================

All outputs go into ``{out}/{sample}/duplex/{ref-1}__{ref-2}/full/``.

``duplex-position-table.csv``
    Per-position data of the fused molecule: the 5' strand's data followed by
    three empty linker positions followed by the 3' strand's data.
    Positions with no data (the linker, a data-less partner, and any position
    outside a strand's table region) are left empty.

``duplex-report.json``
    Which two sources were combined, where the strand break falls, and how
    many clusters the duplex has.
    See :doc:`/formats/report/duplex`.

``refseq.brickle``
    The fused reference sequence, stored like any other reference sequence.


Quick example
================================================================================

Duplex every pair of references that were filtered in one sample::

    seismic duplex out/sample-1/filter

Duplex a target reference with an antisense oligonucleotide, and nothing
else::

    seismic duplex out/sample-1/filter/ref-1 --no-duplex-pair \
        --duplex-sequence ASO ATTCACTTTCATAATGCTGG

Model a homodimer of each input reference::

    seismic duplex out/sample-1/filter --no-duplex-pair --dimer

Then cofold the duplexes, which produces ordinary CT and dot-bracket files::

    seismic fold out/sample-1/duplex


Options
================================================================================

Choosing which duplexes to make
    Each of these contributes its own duplexes, so any combination of them can
    be used at once.

    ``--duplex-pair`` / ``--no-duplex-pair``
        Duplex every pair of two different input tables (default: on).
        Only tables that share the same branches are paired.
    ``--dimer`` / ``--no-dimer``
        Duplex each input table with itself, to model a homodimer
        (default: off).
        A homodimer keeps both names (``ref-1__ref-1``) so that its output is
        never confused with the structure of the reference alone.
    ``--duplex-file FILE``
        Duplex each input table with every sequence in this FASTA file.
        These partners carry no mutational data.
    ``--duplex-sequence NAME SEQUENCE``
        Duplex each input table with one named strand given directly, 5' to
        3', e.g. ``--duplex-sequence ASO ATTCACTTTCATAATGCTGG``.
        Repeat the option to add more partners.  These partners carry no
        mutational data.

Choosing how much of each reference to include
    ``--duplex-full`` / ``--duplex-table-region``
        Whether each strand is the full reference sequence with the table's
        data mapped onto its region (``--duplex-full``, the default), or only
        the part of the reference that the table covers
        (``--duplex-table-region``).
        Use the default when sequence outside the analyzed region can still
        form base pairs; use ``--duplex-table-region`` to restrict folding to
        the region you actually probed.
    ``--duplex-full-ref REF`` and ``--duplex-region-ref REF``
        Override the above for one reference at a time (repeatable).
        Giving the same reference to both raises an error.

Branches
    ``--branch NAME`` (``-b``)
        Write outputs to ``{out}/{sample}/duplex_{NAME}/``.
        See :doc:`/use/branch`.

Performance
    ``--num-cpus N`` — multiprocessing; see :doc:`/use/parallel`.
    ``--force`` — overwrite existing outputs.

The auto-generated :doc:`/cli` lists every option with its current default.


Caveats
================================================================================

- Cofolding a duplex requires RNAcofold_ from ViennaRNA_, regardless of
  ``--fold-backend``; RNAstructure cannot fold duplexes in SEISMIC-RNA.
- Cofolding supports only the ``Deigan`` and ``Cordero`` energy methods
  (``--fold-energy-method``); ``auto`` chooses ``Cordero`` for DMS and
  ``Deigan`` otherwise.
- Each strand's reactivities are normalized separately, because the two
  strands can come from different samples whose overall mutation rates differ.
- ``--duplex-pair`` grows quickly: *n* input tables give *n*(*n* − 1)/2
  duplexes.  Narrow the input paths, or turn it off with
  ``--no-duplex-pair``, when you want only specific combinations.


Common unexpected results
================================================================================

Warning: nothing to duplex
    Only one table was found and no second strand was requested.
    Give at least two tables, or add ``--dimer``, ``--duplex-file``, or
    ``--duplex-sequence``.

Error: cannot duplex references from different branches
    The two strands came from different :doc:`branches </use/branch>`.
    A duplex records one branch chain, so both data-bearing strands must
    share it.  (A partner without data has no branches and always fits.)

Error: cannot duplex references probed with different chemicals
    The duplex is folded with a single energy method, so both data-bearing
    strands must have been probed with the same chemical (set by ``--probe``
    at the :doc:`filter` step).

Warning: duplicate reference, region, and sample
    The same table was given more than once among the inputs; only one copy
    is used.

Warning: no best number of clusters
    A clustered strand whose numbers of clusters all failed the
    :doc:`cluster` step's filters is used at the one number of clusters its
    table contains.  If such a table contains several, which to duplex is
    ambiguous and raises an error; rerun :doc:`cluster` or pick a table with
    a single number of clusters.

More output than expected from clustered inputs
    The duplex of two clustered strands is the cross-product of their
    clusters, so the number of profiles is the product of the two numbers of
    clusters.


See also
================================================================================

- :doc:`filter`, :doc:`cluster` — produce the tables this step combines
- :doc:`fold` — cofolds the duplexes this step makes
- :doc:`draw`, :doc:`graph` — visualize duplex structures and profiles
- :doc:`/formats/report/duplex`
- :doc:`/use/branch`, :doc:`/use/inputs`, :doc:`/use/parallel`


.. _RNAcofold: https://www.tbi.univie.ac.at/RNA/RNAcofold.1.html
.. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/
