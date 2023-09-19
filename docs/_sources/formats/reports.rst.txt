
Report Formats
========================================================================


Relate Report
------------------------------------------------------------------------

Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======================== ========= ===============================================
Name                     Type      Description
======================== ========= ===============================================
Name of Sample           str       Sample's name (directory of XAM file)
Name of Reference        str       Reference sequence's name (name of XAM file)
Sequence of Reference    str       Reference sequence as DNA (from FASTA file)
Length of Sequence (nt)  int       Number of nucleotides in the reference sequence
Number of Reads Passed   int       Number of reads processed into relation vectors
Number of Reads Failed   int       Number of reads that failed to be processed
MD5 Checksums of Batches list[str] List of the MD5 checksum of each parquet file
Number of Batches        int       Number of parquet files
Time Began               str       Date and time at which the step began
Time Ended               str       Date and time at which the step ended
Time Taken (minutes)     float     Duration of the step, in minutes
Speed (reads per minute) float     Speed, in number of reads processed per minute
======================== ========= ===============================================

Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Name of Sample": "sars2-fse",
        "Name of Reference": "sars2",
        "Sequence of Reference": "CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCGTGCGGCACAGGCACTAGTACTGATGTCGTATACAGGGCTTTTGACATCTACAATGATAAAGTAGCTGGTTTTGCTAAATTCCTAAAAACTAATTGTTGTCGCTTCCAAGAAAAGGACGAAG",
        "Length of Sequence (nt)": 283,
        "Number of Reads Passed": 272516,
        "Number of Reads Failed": 0,
        "MD5 Checksums of Batches": [
            "f9874e5c37b56316d47bacc0a7cc604e",
            "15d30bd7b2075d51d7147fc5a454ea1a"
        ],
        "Number of Batches": 2,
        "Time Began": "2023-08-28 at 16:13:11",
        "Time Ended": "2023-08-28 at 16:14:12",
        "Time Taken (minutes)": 1.02,
        "Speed (reads per minute)": 266649.0
    }
