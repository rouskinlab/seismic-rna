
Clusterscan Report
------------------

Clusterscan Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^

========================================================================================================================================= ===============
Name                                                                                                                                      Data Type
========================================================================================================================================= ===============
Sample                                                                                                                                    str
Branches                                                                                                                                  dict[str, str]
Reference                                                                                                                                 str
Region                                                                                                                                    str
Best number of clusters for each domain, keyed by domain coordinates (end5, end3)                                                         dict[str, int]
Filterscan domains (end5, end3) merged into each final domain, keyed by the final domain (end5, end3); omits domains that were not merged dict[str, list]
Time began                                                                                                                                str
Time ended                                                                                                                                str
Time taken (minutes)                                                                                                                      float
Version of SEISMIC-RNA                                                                                                                    str
========================================================================================================================================= ===============

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
        "Best number of clusters for each domain, keyed by domain coordinates (end5, end3)": {
            "3,56": 2,
            "124,230": 3
        },
        "Filterscan domains (end5, end3) merged into each final domain, keyed by the final domain (end5, end3); omits domains that were not merged": {
            "124,230": [
                [124, 174],
                [175, 230]
            ]
        },
        "Time began": "2026-07-02 at 10:07:48",
        "Time ended": "2026-07-02 at 10:14:39",
        "Time taken (minutes)": 6.85,
        "Version of SEISMIC-RNA": "0.26.0dev"
    }
