
Key algorithms
++++++++++++++

Counting mutations for the entire bitvector
-------------------------------------------

The following table describes the aggregation of the bitvector into a set of per-residue counts (except for ``sub_hist`` which is per-read). 

Please refer to :ref:`bitvector` for the definition of the bits and the bitvector format.


======================= =========================================================================================================================
**Aggregation**          **Description**
----------------------- -------------------------------------------------------------------------------------------------------------------------
``cov``                  Per-residue count of covered bases
``del``                  Per-residue count of deleted bases
``ins``                  Per-residue count of inserted bases
``info``                 Per-residue count of bases that are substituted for another base or non-mutated. Doesn't include deleted bases. 
``sub_A``                Per-residue count of bases that substituted to A
``sub_C``                Per-residue count of bases that substituted to C
``sub_G``                Per-residue count of bases that substituted to G
``sub_T``                Per-residue count of bases that substituted to T
``sub_N``                Per-residue count of bases that substituted to any base
``sub_rate``             Per-residue ratio ``sub_N/info``
``sub_hist``             Per-read count of mutations histogram (i.e. number of reads with 0 mutations, 1 mutation, 2 mutations, etc.)
======================= =========================================================================================================================



**Example**:

.. code-block:: bash

    # example.orc
                    # Bitvector 
                    C    A    C    A    A    T    G    T    G   # reference sequence 
                    9    5    1    128  1    0    1    1    2   # read 1
                    1    1    128  1    0    1    16   2    1   # read 2
                    64   1    2    1    1    1    16   1    1   # read 3

        # Read 1 
    cov             1    1    1    1    1    0    1    1    1 
    del             0    0    0    0    0    0    0    0    1
    ins             0    1    0    0    0    0    0    0    0  
    info            1    1    1    1    1    0    1    1    0 
    sub_A           0    0    0    0    0    0    0    0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           0    0    0    0    0    0    0    0    0
    sub_T           0    0    0    1    0    0    0    0    0
    sub_N           0    0    0    1    0    0    0    0    0
    sub_rate        0    0    0    1    0    0    0    0    0

    # Read 2
    cov             1    1    1    1    0    1    1    1    1 
    del             0    0    0    0    0    0    0    1    0
    ins             0    0    0    0    0    0    0    0    0  
    info            1    1    1    1    0    1    1    0    1
    sub_A           0    0    0    0    0    0    1    0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           0    0    0    0    0    0    0    0    0
    sub_T           0    0    1    0    0    0    0    0    0
    sub_N           0    0    1    0    0    0    1    0    0
    sub_rate        0    0    1    0    0    0    1    0    0

    # Read 3
    cov             1    1    1    1    1    1    1    1    1 
    del             0    0    1    0    0    0    0    1    0
    ins             0    0    0    0    0    0    0    0    0  
    info            1    1    0    1    1    1    1    1    1
    sub_A           0    0    0    0    0    0    1    0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           1    0    0    0    0    0    0    0    0
    sub_T           0    0    0    0    0    0    0    0    0
    sub_N           1    0    0    0    0    0    1    0    0
    sub_rate        1    0    0    0    0    0    1    0    0

    # DREEM output = sum of these

    # Aggregation ---------------------------------------------
    cov             3    3    3    3    2    2    3    3    3
    del             0    0    0    0    0    0    0    1    1
    ins             0    1    0    0    0    0    0    0    0  
    info            3    3    2    3    2    2    3    2    2
    sub_A           0    0    0    0    0    0    2    0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           1    0    0    0    0    0    0    0    0
    sub_T           0    0    1    1    0    0    0    0    0
    sub_N           1    0    1    1    0    0    2    0    0
    sub_N           1    0    1    1    0    0    2    0    0
    sub_rate        0.33 0.00 0.5  0.33 0.00 0.00 0.66 0.00 0.00

    # This one is per read, not per base
    sub_hist   0  1  2  0  0  0  0  0  0  0  



Counting mutations for a certain cluster #TODO
----------------------------------------------

Clustering outputs for each read the likelihood to belong to a certain cluster. 
The sum of the likelihoods for different clusters for a read is 1.

When counting mutations for a certain cluster, we weight the mutations by the likelihood of the read to belong to that cluster.

Let's add likelihoods to belong to a cluster K2_1 to the example above:

.. code-block:: bash

    # example.orc
                    # Bitvector 
                    C    A    C    A    A    T    G    T    G   # reference sequence 
                    9    5    1    128  1    0    1    1    2   # read 1, likelihood=1
                    1    1    128  1    0    1    16   2    1   # read 2  likelihood=0.5
                    64   1    2    1    1    1    16   1    1   # read 3  likelihood=0.1

    # Read 1, likelihood=1
    cov             1    1    1    1    1    0    1    1    1 
    del             0    0    0    0    0    0    0    0    1
    ins             0    1    0    0    0    0    0    0    0  
    info            1    1    1    1    1    0    1    1    0 
    sub_A           0    0    0    0    0    0    0    0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           0    0    0    0    0    0    0    0    0
    sub_T           0    0    0    1    0    0    0    0    0
    sub_N           0    0    0    1    0    0    0    0    0
    sub_rate        0    0    0    1    0    0    0    0    0

    # Read 2, likelihood=0.5
    cov             1    1    1    1    0    1    1    1    1 
    del             0    0    0    0    0    0    0    1    0
    ins             0    0    0    0    0    0    0    0    0  
    info            1    1    1    1    0    1    1    0    1
    sub_A           0    0    0    0    0    0    1    0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           0    0    0    0    0    0    0    0    0
    sub_T           0    0    1    0    0    0    0    0    0
    sub_N           0    0    1    0    0    0    1    0    0
    sub_rate        0    0    1    0    0    0    1    0    0

    # Read 3, likelihood=0.1
    cov             1    1    1    1    1    1    1    1    1 
    del             0    0    1    0    0    0    0    0    0
    ins             0    0    0    0    0    0    0    0    0  
    info            1    1    0    1    1    1    1    1    1
    sub_A           0    0    0    0    0    0    1    0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           1    0    0    0    0    0    0    0    0
    sub_T           0    0    0    0    0    0    0    0    0
    sub_N           1    0    0    0    0    0    1    0    0
    sub_rate        1    0    0    0    0    0    1    0    0

    # DREEM output = sum weighted by the likelihood: read1 + 0.5*read2 + 0.1*read3

    # Aggregation ---------------------------------------------
    cov             1.6  1.6  1.6  1.6  1.1  0.6  1.6  1.6  1.6
    del             0    0    0.1  0    0    0    0    0.5  1
    ins             0    1    0    0    0    0    0    0    0  
    info            1.6  1.6  1.5  1.6  1.1  0.6  1.6  1.1  0.6
    sub_A           0    0    0    0    0    0    0.6  0    0
    sub_C           0    0    0    0    0    0    0    0    0
    sub_G           0.1  0    0    0    0    0    0    0    0
    sub_T           0    0    0.5  1    0    0    0    0    0
    sub_N           0.1  0    0.5  1    0    0    0.6  0    0
    sub_rate        0.06 0    0.33 0.62 0    0    0.37 0    0        

    # This one is per read, not per base
    sub_hist  0  1  0.6 

Predicting the structure and the free energy
--------------------------------------------

We use the following script:

.. code-block:: text

    echo ">ref\nGGCGACACAGTCGACGGTTTTCACA">temp.fasta
    Fold temp.fasta temp.ct
    ct2dot temp.ct 1 temp_dot.txt
    cat temp_dot.txt

.. note::

    The results are stored in a temporary json file ``ledger.json``, so that if you predict multiple times the same sequence, you will not have to wait for the prediction to be done again.