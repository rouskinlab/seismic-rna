
Duplex Report
-------------

Duplex Report: Fields
^^^^^^^^^^^^^^^^^^^^^

======================================================== ==============
Name                                                     Data Type     
======================================================== ==============
Sample                                                   str           
Branches                                                 dict[str, str]
Reference                                                str           
Region                                                   str           
Region 5' end                                            int           
Region 3' end                                            int           
Length of the 5' strand (strand-break position)          int           
Numbers of clusters in the duplex (empty if unclustered) list[int]     
Source table files                                       list[str]     
Checksum of reference sequence (SHA-512)                 str           
Time began                                               str           
Time ended                                               str           
Time taken (minutes)                                     float         
Version of SEISMIC-RNA                                   str           
======================================================== ==============

``Reference`` and ``Sample`` name both sources joined by ``__``, and ``Region``
is always ``full`` because a duplex spans its entire fused sequence.
``Length of the 5' strand (strand-break position)`` is where the 5' strand ends;
the three-nucleotide linker follows it, and the 3' strand begins after that.
``Source table files`` are the two source position tables, relative to the
output directory; the entry for a strand given as a bare sequence (from
``--duplex-file`` or ``--duplex-sequence``) is empty, because it has no data.

Duplex Report: Example
^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": "",
            "idmut": "",
            "filter": "",
            "duplex": ""
        },
        "Reference": "refA__refB",
        "Region": "full",
        "Region 5' end": 1,
        "Region 3' end": 123,
        "Length of the 5' strand (strand-break position)": 60,
        "Numbers of clusters in the duplex (empty if unclustered)": [],
        "Source table files": [
            "sample1/filter/refA/full/filter-position-table.csv",
            "sample1/filter/refB/full/filter-position-table.csv"
        ],
        "Checksum of reference sequence (SHA-512)": "ec58bd95183a3cfef8a1200d5d45adef8f51ccf222ae401095f0040d7ca50f4ede3610182d2266268e522e576cd4174c02cc620c21eee5020833d50b368d8466",
        "Time began": "2026-09-02 at 23:25:14",
        "Time ended": "2026-09-02 at 23:25:14",
        "Time taken (minutes)": 0.0,
        "Version of SEISMIC-RNA": "0.26.0"
    }
