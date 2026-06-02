
Pool Report
-----------

Pool Report: Fields
^^^^^^^^^^^^^^^^^^^

===================================================================== ==============
Name                                                                  Data Type     
===================================================================== ==============
Sample                                                                str           
Branches                                                              dict[str, str]
Reference                                                             str           
Pooled samples                                                        list[str]     
Pool samples only if their Pearson correlation is at least this large float         
Pool samples only if their mean arsince distance is at most this      float         
Time began                                                            str           
Time ended                                                            str           
Time taken (minutes)                                                  float         
Version of SEISMIC-RNA                                                str           
===================================================================== ==============

Pool Report: Example
^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "pooled",
        "Branches": {
            "align": "",
            "idmut": ""
        },
        "Reference": "myref",
        "Pooled samples": [
            "sample1",
            "sample2"
        ],
        "Pool samples only if their Pearson correlation is at least this large": 0.0,
        "Pool samples only if their mean arsince distance is at most this": 1.0,
        "Time began": "2026-05-31 at 13:01:12",
        "Time ended": "2026-05-31 at 13:01:15",
        "Time taken (minutes)": 0.04,
        "Version of SEISMIC-RNA": "0.25.3"
    }
