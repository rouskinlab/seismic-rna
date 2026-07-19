
Filterscan Report
-----------------

Filterscan Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^

==================================================================================================================================================================================================================================================================================== ===============
Name                                                                                                                                                                                                                                                                                 Data Type
==================================================================================================================================================================================================================================================================================== ===============
Sample                                                                                                                                                                                                                                                                               str
Branches                                                                                                                                                                                                                                                                             dict[str, str]
Reference                                                                                                                                                                                                                                                                            str
Region                                                                                                                                                                                                                                                                               str
Make each tile this length (if 0, use 2x the median read length)                                                                                                                                                                                                                     int
Make adjacent tiles overlap by at least this fraction of length                                                                                                                                                                                                                      float
Erase the filter reports/batches from the tiling step                                                                                                                                                                                                                                bool
Consider only pairs of positions no more than this far apart when finding domains (if 0, do not restrict pairs beyond what the tile length already bounds)                                                                                                                           int
Analyze only pairs of positions whose number of jointly covering reads is at least this value when finding domains (pairs with less coverage are too noisy to score reliably)                                                                                                        int
Call domains at this false discovery rate (FDR), calibrated by comparing each candidate domain's summed correlation score against chance blocks found in simulated null replicates: higher values call more (and weaker) domains; lower values call fewer, more conservative domains float
Simulate this many per-position independence-null replicates to calibrate the domain-calling penalty at --domain-fdr (at least 1)                                                                                                                                                    int
Seed for the random number generator                                                                                                                                                                                                                                                 int
Cluster only the regions with at least this many positions                                                                                                                                                                                                                           int
If there are gaps between regions to cluster, OMIT (do not cluster) the gaps, INSERT a new region into each gap, or EXPAND the existing regions to fill the gaps                                                                                                                     str
Coordinates of tiles (end5, end3)                                                                                                                                                                                                                                                    list[list[int]]
Number of pairs with chi-square above the null expectation                                                                                                                                                                                                                                   int
Number of domains detected                                                                                                                                                                                                                                                           int
Coordinates of domains (end5, end3)                                                                                                                                                                                                                                                  list[list[int]]
Time began                                                                                                                                                                                                                                                                           str
Time ended                                                                                                                                                                                                                                                                           str
Time taken (minutes)                                                                                                                                                                                                                                                                 float
Version of SEISMIC-RNA                                                                                                                                                                                                                                                               str
==================================================================================================================================================================================================================================================================================== ===============

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
        "Call domains at this false discovery rate (FDR), calibrated by comparing each candidate domain's summed correlation score against chance blocks found in simulated null replicates: higher values call more (and weaker) domains; lower values call fewer, more conservative domains": 0.1,
        "Simulate this many per-position independence-null replicates to calibrate the domain-calling penalty at --domain-fdr (at least 1)": 10,
        "Seed for the random number generator": 0,
        "Cluster only the regions with at least this many positions": 20,
        "If there are gaps between regions to cluster, OMIT (do not cluster) the gaps, INSERT a new region into each gap, or EXPAND the existing regions to fill the gaps": "omit",
        "Coordinates of tiles (end5, end3)": [
            [
                1,
                120
            ]
        ],
        "Number of pairs with chi-square above the null expectation": 270,
        "Number of domains detected": 1,
        "Coordinates of domains (end5, end3)": [
            [
                1,
                56
            ]
        ],
        "Time began": "2026-07-16 at 23:28:48",
        "Time ended": "2026-07-16 at 23:29:19",
        "Time taken (minutes)": 0.52,
        "Version of SEISMIC-RNA": "0.26.0dev"
    }
