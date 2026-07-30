
Clusterscan Report
------------------

Clusterscan Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^

===================================================================== ====================
Name                                                                  Data Type
===================================================================== ====================
Sample                                                                str
Branches                                                              dict[str, str]
Reference                                                             str
Region                                                                str
Coordinates of domains (end5, end3)                                   list[list[int]]
Best number of clusters for each domain                               list[int]
Original filterscan domains (end5, end3) comprising each final domain list
Time began                                                            str
Time ended                                                            str
Time taken (minutes)                                                  float
Version of SEISMIC-RNA                                                str
===================================================================== ====================

Clusterscan Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "sample1",
        "Branches": {
            "idmut": "",
            "clusterscan": ""
        },
        "Reference": "myref",
        "Region": "full",
        "Coordinates of domains (end5, end3)": [
            [3, 56],
            [124, 174]
        ],
        "Best number of clusters for each domain": [
            2,
            2
        ],
        "Original filterscan domains (end5, end3) comprising each final domain": [
            [
                [3, 56]
            ],
            [
                [124, 174]
            ]
        ],
        "Time began": "2026-07-02 at 10:07:48",
        "Time ended": "2026-07-02 at 10:14:39",
        "Time taken (minutes)": 6.85,
        "Version of SEISMIC-RNA": "0.26.0dev"
    }
