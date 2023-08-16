Mutation Vector
+++++++++++++++

Mutation vector files store the mutations in each read as a table in the Apache Optimized Row Columnar (ORC) format (extension ``.orc``).
Each row of the table corresponds to one read and each column to one position in the reference sequence.
Accordingly, each row is labeled on its left side with the name of the read (from the QNAME field in the SAM file), and each column is labeled at the top with the base in the reference sequence (A, C, G, or T) followed by the position in the reference sequence (as in the example below).
The base at the 5' end of the reference sequence is numbered 1, and the base at the 3' end is numbered with the length of the sequence. DREEM allows generating mutation vectors over an entire reference sequence as well as over a section that may begin and end in the middle of the sequence.
The example below shows mutation vectors that begin at position 36 of the reference, and the base at that position is cytosine (C).

**Example**:

An example showing the first three rows (reads) in a mutation vector file. Each element of the matrix is a hexadecimal representation of an integer from 0 (``00``) to 255 (``ff``)::

    Read C36 A37 C38 C39 A40 T41 C42 T43 G44
       A  00  00  01  01  80  01  02  01  01
       B  01  01  03  03  01  01  01  01  01  
       C  01  20  01  01  e1  01  01  00  00


Each element of the matrix contains 8 bits, each of which indicates a possible relationship between the read and the position in the reference sequence, according to the following table:

 ========== ===== ===== ========================================================= 
  Bit        Dec   Hex   Relationship                                                  
 ========== ===== ===== ========================================================= 
  00000001     1    01   match                                                    
  00000010     2    02   deletion                                                 
  00000100     4    04   insertion of ≥1 base(s) immediately 3' of this position  
  00001000     8    08   insertion of ≥1 base(s) immediately 5' of this position  
  00010000    16    10   substitution to A                                        
  00100000    32    20   substitution to C                                        
  01000000    64    40   substitution to G                                        
  10000000   128    80   substitution to T                                        
 ========== ===== ===== ========================================================= 

- Bit: which one of the 8 bits is set to 1
- Dec: decimal (base 10) value of the binary number in the Bit column
- Hex: hexadecimal (base 16) value of the binary number in the Bit column
- Relationship: how the read compares with the position in the reference

Usually, each element of the matrix has at most one bit set to 1 in the binary representation (e.g. ``00010000``).
In the example matrix above, at read B, position 36, the element is ``01`` (hexadecimal) or ``00000001`` (binary), meaning that the base in the read is the same as the base in the reference (which is a C).
In read C at position 37, the element is ``20`` (hexadecimal) or ``00100000`` (binary), meaning that the read has a C at this position (while the reference sequence has an A).

Multiple bits may be set to 1 to indicate where the relationship between the read and the reference is ambiguous. 
For example, in read B, positions 38 and 39 are both marked ``03``, which is ``00000011`` in binary, with both the "match" and "deletion" bits turned on.
This value means that it is ambiguous whether the read contains a match or deletion at this position.
Read C at position 41 has a value of ``e1``, which is ``11100001`` in binary, meaning that it could be a substitution to C, G, or T, or a match (i.e. A).
This situation happens if the base call has a low Phred score so its value cannot be trusted.

If all bits are zero (``00000000``), it means that the read did not align to that position in the reference sequence.
All-zero values occur only on the edges of a read, never in the middle, as can be seen at the first two positions of read A and the last two positions of read C.

See more about the vectorization algorithm :ref:`here <vector_algos>`.
