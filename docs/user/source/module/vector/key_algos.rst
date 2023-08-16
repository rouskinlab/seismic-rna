
Key algorithm(s)
++++++++++++++++

Define batches of reads
-----------------------

Mutation vectors are generated for aligned reads in batches, rather than all at once.
Batching not only reduces the memory footprint but also increases the speed (as batches can be processed in parallel).
The batch size is 32 million bytes (e.g. 100,000 mutation vectors for a 320 nt reference) by default and can be customized.
The aligned reads are batched using the following algorithm:

1. Compute the intended number of vectors per batch, N = B / L (where B is the batch size in bytes and L is the length of reference sequence).
2. Seek to the beginning of the first aligned read in the SAM file.
3. Define the start of the batch as the current position in the file.
4. Read either N lines (files of single-end reads) or 2N - 1 lines (files of paired-end reads).
5. For single-end files: define the end of the batch as the position of the end of the Nth read.
6. For paired-end files: define the end of the batch as the position of the end of the (2N)th read if the (2N - 1)th and (2N)th reads are mated, otherwise the end of the (2N - 1)th read.
7. Output the start and end positions of the batch.
8. Repeat steps 3 - 7 until reaching the end of the SAM file.

Note that the batch will contain more mutation vectors than intended if any paired-end reads have only one mate aligned,
because then the assumption that each paired-end read occupies two lines will be incorrect.
The maximum batch size (when all paired-end reads have only one mate) is twice the intended size.

Parse CIGAR strings into vectors
--------------------------------

A mutation vector is then created for each read in a batch using the following algorithm:

1. Create an array of L bytes, where L is the length of the reference sequence. Initialize every byte to 255 (i.e. 11111111 in binary).
2. Create two variables to track the positions in the read (initialized to 0) and reference (initialized to the left-most position in the reference sequence to which the read aligns).
3. Read one operation from the CIGAR string. For every position covered by the operation:
    a. If the operation is a match (=), then if the read quality is sufficient, mark the position as a match, otherwise an ambiguous match or substitution.
    b. If the operation is a substitution (X) or alignment (M), then if the read quality is sufficient, mark the position as a match if the read and reference agree or a substitution if they disagree, otherwise an ambiguous match or substitution.
    c. If the operation is a deletion (D), then mark the position as a deletion and create a Deletion object (for later).
    d. If the operation is an insertion (I), then create an Insertion object (for later) but do not mark any positions at this time.
    e. If the operation is a soft clipping (S), then skip it.
    f. Otherwise, raise an error. (This behavior may be changed in the future to accommodate more types of CIGAR operations.)
    g. If the operation consumes the read and/or reference, then advance to the next position in the corresponding sequence.
4. For each Insertion object, mark the positions on its immediate 5' and 3' sides as being next to an insertion.

Find ambiguous mutations
------------------------

Insertions and deletions can cause ambiguous alignments.
For example, if the read contains CC where the reference contains CCC, then the read contains a deletion of one C, but which was deleted cannot be determined.
The following algorithm marks all ambiguous positions caused by insertions and deletions:



Merge vectors of paired-end reads
---------------------------------


Write output files
------------------
