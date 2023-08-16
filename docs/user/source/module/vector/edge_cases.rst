
Edge cases handling
+++++++++++++++++++

CIGAR Operations
----------------
For now, vectoring only recognizes the following CIGAR string operations:
matches (=), substitutions (X), alignments (M), insertions (I), deletions (D), and soft clippings (S).
Other operations raise an error.
If needed, future versions of DREEM may handle additional operations that are defined in the SAM format specification, such as skipped reference positions (N) and hard clippings (H).

Ambiguous Insertions and Deletions
----------------------------------
The algorithm that finds ambiguous insertions and deletions assumes that an insertion and a deletion cannot be next to each other
(although two or more insertions may be next to each other, as may two or more deletions).
Handling insertions next to deletions would require a more complex algorithm because an adjacent insertion and deletion could be the equivalent of a substitution.
Thus, the algorithm would need to consider marking ambiguous positions of not only insertions and deletions but also of substitutions.
Instead, the alignment parameters are set so as to make it impossible for an insertion-deletion pair to score better than a substitution, which circumvents this limitation of the algorithm.

Blank Mutation Vectors
----------------------
If, for whatever reason, the mutation vector returned from vectoring is blank, then the vector is discarded.

Zero Mutation Vectors
---------------------
DREEM will handle (without crashing) potential problems caused by outputting zero mutation vectors, such as division-by-zero errors.

Corrupted SAM Files
-------------------
If a single line on a SAM file cannot be parsed into a valid aligned read, then DREEM will log an error and skip the read without crashing.
If, for whatever reason, an entire BAM/SAM file or batch of reads cannot be processed, then DREEM will log a critical error and skip the batch or file.
