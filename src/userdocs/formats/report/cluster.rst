
Cluster Report
--------------

Cluster Report: Fields
^^^^^^^^^^^^^^^^^^^^^^

==================================================================================================== ====================
Name                                                                                                 Data Type           
==================================================================================================== ====================
Sample                                                                                               str                 
Branches                                                                                             dict[str, str]      
Reference                                                                                            str                 
Region                                                                                               str                 
Seed for the random number generator                                                                 int                 
Start at this many clusters                                                                          int                 
Stop at this many clusters (0 for no limit)                                                          int                 
Try all numbers of clusters (Ks), even after finding the best number                                 bool                
Write all numbers of clusters (Ks), rather than only the best number                                 bool                
Remove runs with two clusters more similar than this correlation                                     float               
Remove runs with two clusters that differ by less than this MARCD                                    float               
Remove runs where a cluster differs by more than this ARCD from the ensemble average at any position float               
Remove runs where any cluster's Gini coefficient exceeds this limit                                  float               
Calculate the jackpotting quotient to find over-represented reads                                    bool                
Confidence level for the jackpotting quotient confidence interval                                    float               
Remove runs whose jackpotting quotient exceeds this limit                                            float               
Maximum number of simulations to compute the jackpotting quotient                                    int                 
Skip calculating the jackpotting quotient if the reads × positions exceeds this limit                int                 
Remove Ks with a log likelihood gap larger than this (0 for no limit)                                float               
Remove Ks where every run has less than this correlation vs. the best                                float               
Remove Ks where every run has more than this MARCD vs. the best                                      float               
Run EM (successfully) at least this number of times for each K                                       int                 
Run EM (successfully or not) at most this number of times for each K                                 int                 
Run EM for at least this many iterations                                                             int                 
Run EM for at most this many iterations                                                              int                 
Stop EM when the log likelihood increases by less than this threshold                                float               
Number of unique reads                                                                               int                 
Whether each number of clusters (K) passed filters                                                   dict[str, bool]     
Best number of clusters                                                                              int                 
Numbers of clusters written to batches                                                               list[int]           
Number of batches                                                                                    int                 
Checksums of batches (SHA-512)                                                                       dict[str, list[str]]
Time began                                                                                           str                 
Time ended                                                                                           str                 
Time taken (minutes)                                                                                 float               
Version of SEISMIC-RNA                                                                               str                 
==================================================================================================== ====================

Cluster Report: Example
^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": "",
            "idmut": "",
            "filter": "",
            "cluster": ""
        },
        "Reference": "myref",
        "Region": "full",
        "Seed for the random number generator": 918720378,
        "Start at this many clusters": 1,
        "Stop at this many clusters (0 for no limit)": 2,
        "Try all numbers of clusters (Ks), even after finding the best number": false,
        "Write all numbers of clusters (Ks), rather than only the best number": false,
        "Remove runs with two clusters more similar than this correlation": 0.9,
        "Remove runs with two clusters that differ by less than this MARCD": 0.016,
        "Remove runs where a cluster differs by more than this ARCD from the ensemble average at any position": 0.2,
        "Remove runs where any cluster's Gini coefficient exceeds this limit": 0.667,
        "Calculate the jackpotting quotient to find over-represented reads": true,
        "Confidence level for the jackpotting quotient confidence interval": 0.95,
        "Remove runs whose jackpotting quotient exceeds this limit": 1.1,
        "Maximum number of simulations to compute the jackpotting quotient": 12,
        "Skip calculating the jackpotting quotient if the reads \u00d7 positions exceeds this limit": 268435456,
        "Remove Ks with a log likelihood gap larger than this (0 for no limit)": 0.0,
        "Remove Ks where every run has less than this correlation vs. the best": 0.97,
        "Remove Ks where every run has more than this MARCD vs. the best": 0.008,
        "Run EM (successfully) at least this number of times for each K": 6,
        "Run EM (successfully or not) at most this number of times for each K": 30,
        "Run EM for at least this many iterations": 10,
        "Run EM for at most this many iterations": 500,
        "Stop EM when the log likelihood increases by less than this threshold": 0.37,
        "Number of unique reads": 1713,
        "Whether each number of clusters (K) passed filters": {
            "1": true,
            "2": true
        },
        "Best number of clusters": 1,
        "Numbers of clusters written to batches": [
            1
        ],
        "Number of batches": 1,
        "Checksums of batches (SHA-512)": {
            "cluster": [
                "f18798a185a2b4ab94d96b7a90939a605701d89795fcce7b5a9cd1f5c0106c5f507e91532b82d5129d8b77983200cfbd1c5f97e9d001e07a503393d54f6e4320"
            ]
        },
        "Time began": "2026-05-31 at 12:53:29",
        "Time ended": "2026-05-31 at 12:58:34",
        "Time taken (minutes)": 5.1,
        "Version of SEISMIC-RNA": "0.25.3"
    }
