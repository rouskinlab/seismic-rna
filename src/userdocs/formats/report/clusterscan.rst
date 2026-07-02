
Clusterscan Report
------------------

Clusterscan Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^

======================================= ====================
Name                                    Data Type           
======================================= ====================
Sample                                  str                 
Branches                                dict[str, str]      
Reference                               str                 
Region                                  str                 
Directories of cluster results          list[str]           
Best number of clusters for each domain list[int]           
Time began                              str                 
Time ended                              str                 
Time taken (minutes)                    float               
Version of SEISMIC-RNA                  str                 
======================================= ====================

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
        "Directories of cluster results": [
            "sample1/cluster/myref/3-56",
            "sample1/cluster/myref/124-174"
        ],
        "Best number of clusters for each domain": [
            2,
            2
        ],
        "Time began": "2026-07-02 at 10:07:48",
        "Time ended": "2026-07-02 at 10:14:39",
        "Time taken (minutes)": 6.85,
        "Version of SEISMIC-RNA": "0.26.0dev"
    }
