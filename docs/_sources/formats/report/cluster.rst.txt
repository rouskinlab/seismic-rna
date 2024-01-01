
Cluster Report
--------------------------------------------------------------------------------

Cluster Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======================================== ======================
Name                                     Data Type
======================================== ======================
Name of Sample                           str
Name of Reference                        str
Name of Section                          str
Number of Unique Bit Vectors             int
Maximum Number of Clusters               int
Number of Independent EM Runs            int
Minimum EM Iterations per Cluster        int
Maximum EM Iterations per Cluster        int
Convergence Threshold for Log Likelihood float
Iterations Until Convergence per Run     dict[int, list[int]]
Log Likelihood per Run                   dict[int, list[float]]
Mean Log Likelihood per Order            dict[int, float]
Std. Dev. Log Likelihood per Order       dict[int, float]
Variation of Information per Order       dict[int, float]
Bayesian Information Criterion per Order dict[int, float]
Optimal Number of Clusters               int
Number of Batches                        int
MD5 Checksums of Batches                 dict[str, list[str]]
Branches                                 list[str]
Time Began                               str
Time Ended                               str
Time Taken (minutes)                     float
Version of SEISMIC-RNA                   str
======================================== ======================

Cluster Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Name of Sample": "ldi",
        "Name of Reference": "sars2_1799",
        "Name of Section": "fse",
        "Number of Unique Bit Vectors": 104850,
        "Maximum Number of Clusters": 2,
        "Number of Independent EM Runs": 1,
        "Minimum EM Iterations per Cluster": 10,
        "Maximum EM Iterations per Cluster": 300,
        "Convergence Threshold for Log Likelihood": 0.01,
        "Iterations Until Convergence per Run": {
            "1": [
                2
            ],
            "2": [
                512
            ]
        },
        "Log Likelihood per Run": {
            "1": [
                -3078650.555
            ],
            "2": [
                -3073073.026
            ]
        },
        "Mean Log Likelihood per Order": {
            "1": -3078650.555,
            "2": -3073073.026
        },
        "Std. Dev. Log Likelihood per Order": {
            "1": 0.0,
            "2": 0.0
        },
        "Variation of Information per Order": {
            "1": 0.0,
            "2": 0.0
        },
        "Bayesian Information Criterion per Order": {
            "1": 6158468.699,
            "2": 6148481.23
        },
        "Optimal Number of Clusters": 2,
        "Number of Batches": 16,
        "MD5 Checksums of Batches": {
            "clust": [
                "6a88a23437d10cdb68e29beebf70abb8",
                "2468a801a2686bd198bd726cc90ff30f",
                "3efe0b5aea70954ff37ced0ce509fc40",
                "d4361a74f58bfcaab96d9448eb6d91dd",
                "b8aafa415ba94b6d8d34e3433316a2a3",
                "21b2f8c6cd976858476cc52297e276b3",
                "d120b2d3a8f67458bac23e68aad745ce",
                "4176c9379dfa59a73954469290c90740",
                "4be0122fc2d700f7d75f236cf11dcca3",
                "1323277603f89c04fcb665a34e2790b1",
                "4265f976d74f05d6a9bbacb4e4efe443",
                "c76bc505a886475d93743202a6d61461",
                "36af6fb7b3e8ab41ac569c5f6547b039",
                "778d808c7d34aa3abf0cfcad1c5eb8a6",
                "e007d6017dae45004c11906a6bc073da",
                "c2e57a06bcde2a3dd99f8ed9ce7c4947"
            ]
        },
        "Branches": [],
        "Time Began": "2023-12-18 at 17:52:50",
        "Time Ended": "2023-12-18 at 17:53:33",
        "Time Taken (minutes)": 0.7,
        "Version of SEISMIC-RNA": "0.10.0"
    }
