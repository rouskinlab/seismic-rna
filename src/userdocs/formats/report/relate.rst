
Relate Report
------------------------------------------------------------------------

Relate Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======================================= ====================
Name                                    Data Type
======================================= ====================
Name of Sample                          str
Name of Reference                       str
Number of Reads                         int
Number of Batches                       int
MD5 Checksums of Batches                dict[str, list[str]]
MD5 Checksum of Reference Sequence File str
Branches                                list[str]
Time Began                              str
Time Ended                              str
Time Taken (minutes)                    float
Version of SEISMIC-RNA                  str
======================================= ====================

Relate Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Name of Sample": "yourBigCoolDiscovery",
        "Name of Reference": "yourFavoriteRNA",
        "Number of Reads": 45871,
        "Number of Batches": 2,
        "MD5 Checksums of Batches": {
            "relate": [
                "f9874e5c37b56316d47bacc0a7cc604e",
                "15d30bd7b2075d51d7147fc5a454ea1a"
            ],
            "qnames": [
                "0d8c7e53719d9ba28afab22881d4ce8a",
                "36157068effbcc1ac6d13487b10b2d85",
            ]
        },
        "MD5 Checksum of Reference Sequence File": "42ae5619577059c555c7f7893aa87739",
        "Branches": []
        "Time Began": "2023-08-28 at 16:13:11",
        "Time Ended": "2023-08-28 at 16:14:12",
        "Time Taken (minutes)": 1.02,
        "Version of SEISMIC-RNA": "0.9.3"
    }
