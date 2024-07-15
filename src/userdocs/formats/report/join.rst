
Join Report
--------------------------------------------------------------------------------

Join Mask Report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Join Mask Report: Fields
================================================================================

====================== =========
Name                   Data Type
====================== =========
Name of Sample         str
Name of Reference      str
Name of Section        str
Joined Sections        list[str]
Branches               list[str]
Time Began             str
Time Ended             str
Time Taken (minutes)   float
Version of SEISMIC-RNA str
====================== =========

Join Mask Report: Example
================================================================================

::

    {
        "Name of Sample": "ldi",
        "Name of Reference": "sars2_1799",
        "Name of Section": "joined",
        "Joined Sections": [
            "partner",
            "thefse"
        ],
        "Branches": [],
        "Time Began": "2024-01-18 at 13:47:51",
        "Time Ended": "2024-01-18 at 13:47:51",
        "Time Taken (minutes)": 0.0,
        "Version of SEISMIC-RNA": "0.13.0"
    }

Join Cluster Report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Join Cluster Report: Fields
================================================================================

====================== ====================================
Name                   Data Type
====================== ====================================
Name of Sample         str
Name of Reference      str
Name of Section        str
Joined Sections        list[str]
Joined Clusters        dict[str, dict[int, dict[int, int]]]
Branches               list[str]
Time Began             str
Time Ended             str
Time Taken (minutes)   float
Version of SEISMIC-RNA str
====================== ====================================

Join Cluster Report: Example
================================================================================

::

    {
        "Name of Sample": "ldi",
        "Name of Reference": "sars2_1799",
        "Name of Section": "joined",
        "Joined Sections": [
            "partner",
            "thefse"
        ],
        "Joined Clusters": {
            "mysection": {
                "1": {
                    "1": 1
                },
                "2": {
                    "1": 1,
                    "2": 2
                },
                "3": {
                    "1": 1,
                    "2": 3,
                    "3": 2
                }
            },
            "yoursection": {
                "1": {
                    "1": 1
                },
                "2": {
                    "1": 2,
                    "2": 1
                },
                "3": {
                    "1": 3,
                    "2": 2,
                    "3": 1
                }
            }
        },
        "Branches": [],
        "Time Began": "2024-01-18 at 16:01:25",
        "Time Ended": "2024-01-18 at 16:01:25",
        "Time Taken (minutes)": 0.0,
        "Version of SEISMIC-RNA": "0.13.0"
    }
