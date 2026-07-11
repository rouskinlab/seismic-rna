
Filterscan Report
-----------------

Filterscan Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^

=========================================================================================================================================================================================================================================================================== ===============
Name                                                                                                                                                                                                                                                                        Data Type
=========================================================================================================================================================================================================================================================================== ===============
Sample                                                                                                                                                                                                                                                                      str
Branches                                                                                                                                                                                                                                                                    dict[str, str]
Reference                                                                                                                                                                                                                                                                   str
Region                                                                                                                                                                                                                                                                      str
Make each tile this length (if 0, use 2x the median read length)                                                                                                                                                                                                            int
Make adjacent tiles overlap by at least this fraction of length                                                                                                                                                                                                             float
Erase the filter reports/batches from the tiling step                                                                                                                                                                                                                       bool
Find correlated pairs at this false discovery rate (FDR)                                                                                                                                                                                                                    float
Cluster only the regions with at least this many correlated pairs                                                                                                                                                                                                           int
Among pairs that survive the endpoint-peak filter, drop any pair whose L1 (Manhattan) distance to its nearest surviving neighbor exceeds this percentile of all such distances. Pairs more isolated than this threshold are treated as noise.                               float
Minimum number of other surviving pairs that must lie within the pair-distance-percentile L1 threshold for a pair to be kept. Setting this above 1 filters out small coincidental clusters of noise pairs ('buddy noise') at the cost of potentially clipping domain edges. int
Cluster only the regions with at least this many positions                                                                                                                                                                                                                  int
If there are gaps between regions to cluster, OMIT (do not cluster) the gaps, INSERT a new region into each gap, or EXPAND the existing regions to fill the gaps                                                                                                            str
Coordinates of tiles (end5, end3)                                                                                                                                                                                                                                           list[list[int]]
Number of significant pairs detected                                                                                                                                                                                                                                        int
Number of domains detected                                                                                                                                                                                                                                                  int
Coordinates of domains (end5, end3)                                                                                                                                                                                                                                         list[list[int]]
Time began                                                                                                                                                                                                                                                                  str
Time ended                                                                                                                                                                                                                                                                  str
Time taken (minutes)                                                                                                                                                                                                                                                        float
Version of SEISMIC-RNA                                                                                                                                                                                                                                                      str
=========================================================================================================================================================================================================================================================================== ===============

Filterscan Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "sample1",
        "Branches": {
            "idmut": "",
            "filterscan": ""
        },
        "Reference": "myref",
        "Region": "full",
        "Make each tile this length (if 0, use 2x the median read length)": 0,
        "Make adjacent tiles overlap by at least this fraction of length": 0.5,
        "Erase the filter reports/batches from the tiling step": true,
        "Find correlated pairs at this false discovery rate (FDR)": 0.05,
        "Cluster only the regions with at least this many correlated pairs": 2,
        "Among pairs that survive the endpoint-peak filter, drop any pair whose L1 (Manhattan) distance to its nearest surviving neighbor exceeds this percentile of all such distances. Pairs more isolated than this threshold are treated as noise.": 95.0,
        "Minimum number of other surviving pairs that must lie within the pair-distance-percentile L1 threshold for a pair to be kept. Setting this above 1 filters out small coincidental clusters of noise pairs ('buddy noise') at the cost of potentially clipping domain edges.": 2,
        "Cluster only the regions with at least this many positions": 20,
        "If there are gaps between regions to cluster, OMIT (do not cluster) the gaps, INSERT a new region into each gap, or EXPAND the existing regions to fill the gaps": "omit",
        "Coordinates of tiles (end5, end3)": [
            [
                1,
                180
            ]
        ],
        "Number of significant pairs detected": 363,
        "Number of domains detected": 2,
        "Coordinates of domains (end5, end3)": [
            [
                3,
                56
            ],
            [
                124,
                174
            ]
        ],
        "Time began": "2026-07-02 at 10:07:22",
        "Time ended": "2026-07-02 at 10:07:48",
        "Time taken (minutes)": 0.44,
        "Version of SEISMIC-RNA": "0.26.0dev"
    }
