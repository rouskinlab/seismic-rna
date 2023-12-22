
Mask Report
--------------------------------------------------------------------------------

Mask Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

====================================================== ====================
Name                                                   Data Type
====================================================== ====================
Name of Sample                                         str
Name of Reference                                      str
Name of Section                                        str
5' end of Section                                      int
3' end of Section                                      int
Count the Following as Mutations                       list[str]
Count the Following as Matches                         list[str]
Exclude G/U Bases                                      bool
Exclude Poly(A) Sequences of at Least This Length (nt) int
Exclude User-Defined Positions                         list[int]
Minimum Number of Informative Reads per Position       int
Maximum Fraction of Mutations per Position             float
Number of Positions Initially Given                    int
Number of Positions Cut -- G/U Base                    int
Number of Positions Cut -- Poly(A) Sequence            int
Number of Positions Cut -- User-Specified              int
Number of Positions Cut -- Too Few Informative Reads   int
Number of Positions Cut -- Too Many Mutations          int
Number of Positions Ultimately Kept                    int
Positions Cut -- G/U Base                              list[int]
Positions Cut -- Poly(A) Sequence                      list[int]
Positions Cut -- User-Specified                        list[int]
Positions Cut -- Too Few Informative Reads             list[int]
Positions Cut -- Too Many Mutations                    list[int]
Positions Ultimately Kept                              list[int]
Minimum Fraction of Informative Positions per Read     float
Maximum Fraction of Mutations per Read                 float
Minimum Gap Between Mutations (nt)                     int
Number of Reads Initially Given                        int
Number of Reads Cut -- Too Few Informative Positions   int
Number of Reads Cut -- Too Many Mutations              int
Number of Reads Cut -- Mutations Too Close Together    int
Number of Reads Ultimately Kept                        int
Number of Batches                                      int
MD5 Checksums of Batches                               dict[str, list[str]]
Branches                                               list[str]
Time Began                                             str
Time Ended                                             str
Time Taken (minutes)                                   float
Version of SEISMIC-RNA                                 str
====================================================== ====================

Mask Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Name of Sample": "ldi",
        "Name of Reference": "sars2_1799",
        "Name of Section": "fse",
        "5' end of Section": 397,
        "3' end of Section": 426,
        "Count the Following as Mutations": [
            "A -> C",
            "A -> G",
            "A -> T",
            "C -> A",
            "C -> G",
            "C -> T",
            "G -> A",
            "G -> C",
            "G -> T",
            "T -> A",
            "T -> C",
            "T -> G"
        ],
        "Count the Following as Matches": [
            "A -> A",
            "C -> C",
            "G -> G",
            "T -> T"
        ],
        "Exclude G/U Bases": true,
        "Exclude Poly(A) Sequences of at Least This Length (nt)": 5,
        "Exclude User-Defined Positions": [],
        "Minimum Number of Informative Reads per Position": 1000,
        "Maximum Fraction of Mutations per Position": 0.5,
        "Number of Positions Initially Given": 30,
        "Number of Positions Cut -- G/U Base": 15,
        "Number of Positions Cut -- Poly(A) Sequence": 5,
        "Number of Positions Cut -- User-Specified": 0,
        "Number of Positions Cut -- Too Few Informative Reads": 0,
        "Number of Positions Cut -- Too Many Mutations": 0,
        "Number of Positions Ultimately Kept": 10,
        "Positions Cut -- G/U Base": [
            398,
            399,
            400,
            401,
            402,
            403,
            404,
            405,
            407,
            411,
            412,
            415,
            422,
            425,
            426
        ],
        "Positions Cut -- Poly(A) Sequence": [
            416,
            417,
            418,
            419,
            420
        ],
        "Positions Cut -- User-Specified": [],
        "Positions Cut -- Too Few Informative Reads": [],
        "Positions Cut -- Too Many Mutations": [],
        "Positions Ultimately Kept": [
            397,
            406,
            408,
            409,
            410,
            413,
            414,
            421,
            423,
            424
        ],
        "Minimum Fraction of Informative Positions per Read": 0.0,
        "Maximum Fraction of Mutations per Read": 0.1,
        "Minimum Gap Between Mutations (nt)": 0,
        "Number of Reads Initially Given": 562325,
        "Number of Reads Cut -- Too Few Informative Positions": 0,
        "Number of Reads Cut -- Too Many Mutations": 289974,
        "Number of Reads Cut -- Mutations Too Close Together": 0,
        "Number of Reads Ultimately Kept": 272351,
        "Number of Batches": 16,
        "MD5 Checksums of Batches": {
            "mask": [
                "d295b42cc757255ad29070f00fb6a51d",
                "9f8d6ba42947fb1ee3b9a0b67879e99e",
                "194824b402bced072049a46d8c7effb6",
                "b723dc6f19f93bca0827155cca36965d",
                "dfc527236af3a5888d3753b2c0041dc3",
                "27eb008e907d48d932626098f4587178",
                "c2c0ee505c970d55d2e17e8a2d8ba48c",
                "42babd92a5ffdd13e29441f034a6b167",
                "901b083617e6a66df5a68dc03323b791",
                "a13183811fadb0a9906d71a729aee04d",
                "bc18c08eca71f48d57e0119392042cc3",
                "1680be76740c49f0c216574ca61f3ed4",
                "beebe1e830d72da7296aeb62a29e057b",
                "81f2ed98056df9611354e0aa3a219971",
                "09c4d2b19c4c720bd24b2f432b908cb9",
                "c5cec619b953ff80609f172a04abf6ec"
            ]
        },
        "Branches": [],
        "Time Began": "2023-12-18 at 17:51:53",
        "Time Ended": "2023-12-18 at 17:52:09",
        "Time Taken (minutes)": 0.26,
        "Version of SEISMIC-RNA": "0.10.0"
    }
