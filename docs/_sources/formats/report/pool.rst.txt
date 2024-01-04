
Pool Report
------------------------------------------------------------------------

Pool Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

====================== =========
Name                   Data Type
====================== =========
Name of Sample         str
Name of Reference      str
Pooled Samples         list[str]
Branches               list[str]
Time Began             str
Time Ended             str
Time Taken (minutes)   float
Version of SEISMIC-RNA str
====================== =========

Pool Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Name of Sample": "mypool",
        "Name of Reference": "sars2_1799",
        "Pooled Samples": [
            "sample1",
            "sample2"
        ],
        "Branches": [],
        "Time Began": "2024-01-04 at 13:39:05",
        "Time Ended": "2024-01-04 at 13:39:05",
        "Time Taken (minutes)": 0.0,
        "Version of SEISMIC-RNA": "0.11.0"
    }
