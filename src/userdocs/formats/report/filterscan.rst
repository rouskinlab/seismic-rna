
Filterscan Report
-----------------

Filterscan Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^

================================================================================================================================================================================================================================================================================================================== ===============
Name                                                                                                                                                                                                                                                                                                               Data Type
================================================================================================================================================================================================================================================================================================================== ===============
Sample                                                                                                                                                                                                                                                                                                             str
Branches                                                                                                                                                                                                                                                                                                           dict[str, str]
Reference                                                                                                                                                                                                                                                                                                          str
Region                                                                                                                                                                                                                                                                                                             str
Make each tile this length (if 0, use 2x the median read length)                                                                                                                                                                                                                                                   int
Make adjacent tiles overlap by at least this fraction of length                                                                                                                                                                                                                                                    float
Erase the filter reports/batches from the tiling step                                                                                                                                                                                                                                                              bool
Consider only pairs of positions no more than this far apart when finding domains (if 0, do not restrict pairs beyond what the tile length already bounds)                                                                                                                                                         int
Analyze only pairs of positions whose number of jointly covering reads is at least this value when finding domains (pairs with less coverage are too noisy to score reliably)                                                                                                                                      int
Analyze only pairs of positions whose expected number of reads mutated at both positions (under the assumption that the positions mutate independently) is at least this value when finding domains (standard practice for a chi-square test, which becomes unreliable when an expected count drops below about 5) int
Detect domains at this false discovery rate (FDR), Benjamini-Hochberg-adjusted over every candidate block's exact chi-square p-value: higher values call more (and weaker) domains; lower values call fewer, more conservative domains                                                                             float
Merge adjacent domains separated by a gap whose crossing pairs clear this false discovery rate (FDR), Benjamini-Hochberg-adjusted over every cut in the gap: higher values merge more (and weaker) connections; lower values merge fewer, more conservative ones                                                   float
Keep only the domains with at least this many positions                                                                                                                                                                                                                                                            int
If there are gaps between regions to cluster, OMIT (do not cluster) the gaps, INSERT a new region into each gap, or EXPAND the existing regions to fill the gaps                                                                                                                                                   str
Coordinates of tiles (end5, end3)                                                                                                                                                                                                                                                                                  list[list[int]]
Number of pairs with chi-square above the null expectation                                                                                                                                                                                                                                                         int
Number of domains detected                                                                                                                                                                                                                                                                                         int
Coordinates of domains (end5, end3)                                                                                                                                                                                                                                                                                list[list[int]]
Time began                                                                                                                                                                                                                                                                                                         str
Time ended                                                                                                                                                                                                                                                                                                         str
Time taken (minutes)                                                                                                                                                                                                                                                                                               float
Version of SEISMIC-RNA                                                                                                                                                                                                                                                                                             str
================================================================================================================================================================================================================================================================================================================== ===============

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
        "Consider only pairs of positions no more than this far apart when finding domains (if 0, do not restrict pairs beyond what the tile length already bounds)": 0,
        "Analyze only pairs of positions whose number of jointly covering reads is at least this value when finding domains (pairs with less coverage are too noisy to score reliably)": 1000,
        "Analyze only pairs of positions whose expected number of reads mutated at both positions (under the assumption that the positions mutate independently) is at least this value when finding domains (standard practice for a chi-square test, which becomes unreliable when an expected count drops below about 5)": 5,
        "Detect domains at this false discovery rate (FDR), Benjamini-Hochberg-adjusted over every candidate block's exact chi-square p-value: higher values call more (and weaker) domains; lower values call fewer, more conservative domains": 0.1,
        "Merge adjacent domains separated by a gap whose crossing pairs clear this false discovery rate (FDR), Benjamini-Hochberg-adjusted over every cut in the gap: higher values merge more (and weaker) connections; lower values merge fewer, more conservative ones": 0.1,
        "Keep only the domains with at least this many positions": 20,
        "If there are gaps between regions to cluster, OMIT (do not cluster) the gaps, INSERT a new region into each gap, or EXPAND the existing regions to fill the gaps": "omit",
        "Coordinates of tiles (end5, end3)": [
            [
                1,
                180
            ]
        ],
        "Number of pairs with chi-square above the null expectation": 1301,
        "Number of domains detected": 2,
        "Coordinates of domains (end5, end3)": [
            [
                1,
                56
            ],
            [
                125,
                172
            ]
        ],
        "Time began": "2026-07-18 at 23:49:18",
        "Time ended": "2026-07-18 at 23:50:11",
        "Time taken (minutes)": 0.88,
        "Version of SEISMIC-RNA": "0.26.0dev"
    }
