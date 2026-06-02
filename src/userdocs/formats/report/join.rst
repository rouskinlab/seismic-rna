
Join Report
-----------

Join Filter Report
^^^^^^^^^^^^^^^^^^

Join Filter Report: Fields
==========================

====================== ==============
Name                   Data Type     
====================== ==============
Sample                 str           
Branches               dict[str, str]
Reference              str           
Region                 str           
Joined regions         list[str]     
Time began             str           
Time ended             str           
Time taken (minutes)   float         
Version of SEISMIC-RNA str           
====================== ==============

Join Filter Report: Example
===========================

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": "",
            "idmut": "",
            "filter": ""
        },
        "Reference": "myref",
        "Region": "joinedreg",
        "Joined regions": [
            "1-40",
            "41-80"
        ],
        "Time began": "2026-05-31 at 13:01:38",
        "Time ended": "2026-05-31 at 13:01:51",
        "Time taken (minutes)": 0.23,
        "Version of SEISMIC-RNA": "0.25.3"
    }

Join Cluster Report
^^^^^^^^^^^^^^^^^^^

Join Cluster Report: Fields
===========================

====================== ====================================
Name                   Data Type                           
====================== ====================================
Sample                 str                                 
Branches               dict[str, str]                      
Reference              str                                 
Region                 str                                 
Joined regions         list[str]                           
Joined clusters        dict[str, dict[str, dict[str, int]]]
Time began             str                                 
Time ended             str                                 
Time taken (minutes)   float                               
Version of SEISMIC-RNA str                                 
====================== ====================================

Join Cluster Report: Example
============================

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
        "Region": "joinedclust",
        "Joined regions": [
            "1-40",
            "41-80"
        ],
        "Joined clusters": {
            "1-40": {
                "1": {
                    "1": 1
                }
            },
            "41-80": {
                "1": {
                    "1": 1
                }
            }
        },
        "Time began": "2026-05-31 at 13:08:33",
        "Time ended": "2026-05-31 at 13:08:35",
        "Time taken (minutes)": 0.03,
        "Version of SEISMIC-RNA": "0.25.3"
    }
