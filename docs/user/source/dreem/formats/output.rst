

DREEM output
++++++++++++

DREEM outputs one file per sample. 
The file is a JSON file with the following structure.
See :ref:`aggregate_key_algos` for an explicit description of the cluster-level attributes. 

.. note::

   The same data is used in the previous examples.

::


    {
        # Attributes for the entire sample
        "sample": "01_1_S22_reads", # sample name
        "user": "Lauren",      # who did the experiment
        "date": "1/5/23",     # date of experiment
        "exp_env": "in_vitro", # experiment environment (in_vitro or in_vivo)
        "temperature_k": 310,  # incubation temperature in Kelvin 
        "inc_time_tot_secs": 300, # incubation time in seconds
        "DMS_conc_mM": 10.5,  # DMS concentration in mM
        "buffer": "buffer",   # buffer name for in vitro experiments
        "cell_line": NaN,     # cell line name for in vivo experiments

        # Attributes for each reference
        "3042-O-flank_1=hp1-DB": # reference name
        { 
            "num_aligned": 10147, # number of aligned reads
            "sequence": "TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA",
            "barcode_start": 140, # barcode start position (if you use a barcode) (1-based, inclusive)
            "barcode_end": 151,   # barcode end position (if you use a barcode) (1-based, inclusive)
            "family": "hp1",      # some random attribute for the reference
            
            # Attributes for each section
            "section1": # section name
            {
                "section_start": 69, # section start position (1-based, inclusive) 
                "section_end": 86,   # section end position (1-based, inclusive)
                "deltaG": -9.8, # deltaG for the section (Gibbs free energy)
                "sequence": "CACAGTCGAAAGACTGTG", # sequence for the section
                "structure": "(((((((....)))))))",  # structure for the section, predicted suing RNAstructure
                
                # Attributes for each cluster (see Modules/Aggregate/key algos/counting mutations).
                "pop_avg": # Cluster name. "pop_avg" is the population average, i.e. the un-clustered bitvector.  
                           # Clusters are named K(total # of clusters)_(cluster #), e.g. K2_1, K2_2, K3_1, K3_2, K3_3, etc.
                {
                    # Per-residue count of covered bases
                    "cov": [10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147, 10147], 
                    "min_cov": 10147,
                    # Per-residue count of deleted bases
                    "del": [5, 135, 17, 16, 34, 232, 43, 26, 1, 4, 4, 11, 164, 83, 48, 68, 50, 72], 
                    # Per-residue count of inserted bases
                    "ins": [2, 4, 14, 6, 15, 19, 11, 15, 18, 50, 38, 19, 9, 8, 6, 3, 13, 1], 
                    # Per-residue count of bases that are substituted for another base or non-mutated. Doesn't include deleted bases. 
                    "info": [10003, 9904, 10015, 9969, 9957, 9708, 9979, 9974, 10013, 9935, 9978, 10025, 9639, 9997, 10019, 9937, 9999, 9959],
                    # Per-residue count of bases that substituted to A
                    "sub_A": [23, 0, 11, 0, 41, 3, 14, 61, 0, 0, 0, 71, 0, 5, 4, 56, 4, 62], 
                    # Per-residue count of bases that substituted to C
                    "sub_C": [0, 4, 0, 5, 1, 8, 0, 3, 12, 19, 14, 1, 3, 0, 6, 3, 8, 1],
                    # Per-residue count of bases that substituted to G
                    "sub_G": [10, 2, 5, 0, 0, 2, 3, 0, 79, 4, 5, 0, 11, 6, 4, 0, 0, 0],
                    # Per-residue count of bases that substituted to T
                    "sub_T": [15, 19, 11, 9, 7, 0, 4, 13, 10, 64, 66, 21, 7, 5, 0, 16, 0, 15],
                    # Per-residue count of bases that substituted to any base
                    "sub_N": [48, 26, 27, 14, 49, 13, 21, 77, 101, 87, 85, 93, 21, 16, 14, 75, 12, 78], 
                    # sub_N / info for each residue
                    "sub_rate": [0.004798560431870439, 0.0026252019386106625, 0.0026959560659011485, 0.001404353495837095, 0.004921160992266747, 0.0013391017717346518, 0.002104419280489027, 0.007720072187687989, 0.01008688704683911, 0.00875691997986915, 0.008518741230707557, 0.009276807980049876, 0.002178649237472767, 0.001600480144043213, 0.0013973450444156104, 0.007547549562242125, 0.0012001200120012002, 0.007832111657796967], 
                    # Per-read count of mutations histogram (i.e. number of reads with 0 mutations, 1 mutation, 2 mutations, etc.)
                    "sub_hist": [8199, 1415, 301, 115, 55, 27, 17, 4, 6, 1, 4, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    }
                }
            }

            # Attributes for another section
            "section2":
            {
                ...
            }
        }