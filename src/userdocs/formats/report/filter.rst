
Filter Report
-------------

Filter Report: Fields
^^^^^^^^^^^^^^^^^^^^^

================================================================================================================================ ====================
Name                                                                                                                             Data Type           
================================================================================================================================ ====================
Sample                                                                                                                           str                 
Branches                                                                                                                         dict[str, str]      
Reference                                                                                                                        str                 
Region                                                                                                                           str                 
Region 5' end                                                                                                                    int                 
Region 3' end                                                                                                                    int                 
Use the default options for this chemical probe                                                                                  str                 
Count as mutations                                                                                                               list[str]           
Count as matches                                                                                                                 list[str]           
Mask positions with base A                                                                                                       bool                
Mask positions with base C                                                                                                       bool                
Mask positions with base G                                                                                                       bool                
Mask positions with base U                                                                                                       bool                
Mask stretches of at least this many consecutive A bases (0 disables)                                                            int                 
Mask additional positions from a list                                                                                            list[int]           
Mask positions with fewer than this many informative base calls                                                                  int                 
Mask positions with more than this fraction of mutated base calls                                                                float               
Drop reads with fewer than this many bases covering the region                                                                   int                 
Drop reads covering less than this fraction of the region                                                                        float               
Drop paired-end reads with discontiguous mates                                                                                   bool                
Drop reads with less than this fraction of informative base calls                                                                float               
Drop reads with more than this fraction of mutated base calls                                                                    float               
Filter out mutations separated by fewer than this many bases                                                                     int                 
If two mutations are closer than --min-mut-gap positions, MERGE the mutations, DROP the read, or AUTO-select based on the probe. str                 
Stop the filter step after this many iterations (0 for no limit)                                                                 int                 
Correct observer bias using a quick (typically linear time) heuristic                                                            bool                
Treat mutated fractions under this threshold as 0 with --quick-unbias                                                            float               
Total number of positions in the region                                                                                          int                 
Number of positions masked for having base A                                                                                     int                 
Number of positions masked for having base C                                                                                     int                 
Number of positions masked for having base G                                                                                     int                 
Number of positions masked for having base U                                                                                     int                 
Number of positions masked for having base N                                                                                     int                 
Number of positions in stretches of consecutive A bases                                                                          int                 
Number of positions masked from a list                                                                                           int                 
Number of positions with too few informative base calls                                                                          int                 
Number of positions with too many mutations                                                                                      int                 
Number of positions kept after filtering                                                                                         int                 
Positions masked for having base A                                                                                               list[int]           
Positions masked for having base C                                                                                               list[int]           
Positions masked for having base G                                                                                               list[int]           
Positions masked for having base U                                                                                               list[int]           
Positions masked for having base N                                                                                               list[int]           
Positions in stretches of consecutive A bases                                                                                    list[int]           
Positions masked from a list                                                                                                     list[int]           
Positions with too few informative base calls                                                                                    list[int]           
Positions with too many mutations                                                                                                list[int]           
Positions kept after filtering                                                                                                   list[int]           
Total number of reads before filtering                                                                                           int                 
Number of reads dropped from a list                                                                                              int                 
Number of reads dropped due to too few bases covering the region                                                                 int                 
Number of reads dropped due to covering too small a fraction of the region                                                       int                 
Number of reads dropped due to discontiguous mates                                                                               int                 
Number of reads dropped due to too few informative base calls                                                                    int                 
Number of reads dropped due to too many mutations                                                                                int                 
Number of reads dropped due to two mutations too close                                                                           int                 
Number of reads kept after filtering                                                                                             int                 
Number of iterations until convergence (0 if not converged)                                                                      int                 
Number of batches                                                                                                                int                 
Checksums of batches (SHA-512)                                                                                                   dict[str, list[str]]
Time began                                                                                                                       str                 
Time ended                                                                                                                       str                 
Time taken (minutes)                                                                                                             float               
Version of SEISMIC-RNA                                                                                                           str                 
================================================================================================================================ ====================

Filter Report: Example
^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": "",
            "idmut": "",
            "filter": ""
        },
        "Reference": "myref",
        "Region": "full",
        "Region 5' end": 1,
        "Region 3' end": 80,
        "Use the default options for this chemical probe": "DMS",
        "Count as mutations": [
            "A -> C",
            "A -> D",
            "A -> G",
            "A -> I",
            "A -> T",
            "C -> A",
            "C -> D",
            "C -> G",
            "C -> I",
            "C -> T",
            "G -> A",
            "G -> C",
            "G -> D",
            "G -> I",
            "G -> T",
            "N -> D",
            "N -> I",
            "T -> A",
            "T -> C",
            "T -> D",
            "T -> G",
            "T -> I"
        ],
        "Count as matches": [
            "A -> A",
            "C -> C",
            "G -> G",
            "T -> T"
        ],
        "Mask positions with base A": false,
        "Mask positions with base C": false,
        "Mask positions with base G": true,
        "Mask positions with base U": true,
        "Mask stretches of at least this many consecutive A bases (0 disables)": 5,
        "Mask additional positions from a list": [],
        "Mask positions with fewer than this many informative base calls": 1000,
        "Mask positions with more than this fraction of mutated base calls": 1.0,
        "Drop reads with fewer than this many bases covering the region": 1,
        "Drop reads covering less than this fraction of the region": 0.0,
        "Drop paired-end reads with discontiguous mates": true,
        "Drop reads with less than this fraction of informative base calls": 0.95,
        "Drop reads with more than this fraction of mutated base calls": 1.0,
        "Filter out mutations separated by fewer than this many bases": 4,
        "If two mutations are closer than --min-mut-gap positions, MERGE the mutations, DROP the read, or AUTO-select based on the probe.": "drop",
        "Stop the filter step after this many iterations (0 for no limit)": 0,
        "Correct observer bias using a quick (typically linear time) heuristic": true,
        "Treat mutated fractions under this threshold as 0 with --quick-unbias": 0.001,
        "Total number of positions in the region": 80,
        "Number of positions masked for having base A": 0,
        "Number of positions masked for having base C": 0,
        "Number of positions masked for having base G": 13,
        "Number of positions masked for having base U": 22,
        "Number of positions masked for having base N": 0,
        "Number of positions in stretches of consecutive A bases": 5,
        "Number of positions masked from a list": 0,
        "Number of positions with too few informative base calls": 14,
        "Number of positions with too many mutations": 0,
        "Number of positions kept after filtering": 26,
        "Positions masked for having base A": [],
        "Positions masked for having base C": [],
        "Positions masked for having base G": [
            1,
            6,
            7,
            14,
            31,
            33,
            34,
            35,
            39,
            47,
            49,
            60,
            69
        ],
        "Positions masked for having base U": [
            13,
            15,
            18,
            22,
            25,
            26,
            27,
            30,
            41,
            46,
            51,
            52,
            53,
            54,
            58,
            61,
            62,
            64,
            67,
            71,
            76,
            79
        ],
        "Positions masked for having base N": [],
        "Positions in stretches of consecutive A bases": [
            8,
            9,
            10,
            11,
            12
        ],
        "Positions masked from a list": [],
        "Positions with too few informative base calls": [
            2,
            3,
            4,
            5,
            66,
            68,
            70,
            72,
            73,
            74,
            75,
            77,
            78,
            80
        ],
        "Positions with too many mutations": [],
        "Positions kept after filtering": [
            16,
            17,
            19,
            20,
            21,
            23,
            24,
            28,
            29,
            32,
            36,
            37,
            38,
            40,
            42,
            43,
            44,
            45,
            48,
            50,
            55,
            56,
            57,
            59,
            63,
            65
        ],
        "Total number of reads before filtering": 3595,
        "Number of reads dropped from a list": 0,
        "Number of reads dropped due to too few bases covering the region": 0,
        "Number of reads dropped due to covering too small a fraction of the region": 0,
        "Number of reads dropped due to discontiguous mates": 0,
        "Number of reads dropped due to too few informative base calls": 134,
        "Number of reads dropped due to too many mutations": 0,
        "Number of reads dropped due to two mutations too close": 16,
        "Number of reads kept after filtering": 3445,
        "Number of iterations until convergence (0 if not converged)": 2,
        "Number of batches": 1,
        "Checksums of batches (SHA-512)": {
            "filter": [
                "bcfd01c9340c63d1474e58e72657e845a6a0681988a68aaf205b30d2003e781ee492327f1ea6f44637853309ce323d4714b0458ed3915ddca1ec951751d4c3e6"
            ]
        },
        "Time began": "2026-05-31 at 12:53:14",
        "Time ended": "2026-05-31 at 12:53:26",
        "Time taken (minutes)": 0.2,
        "Version of SEISMIC-RNA": "0.25.3"
    }
