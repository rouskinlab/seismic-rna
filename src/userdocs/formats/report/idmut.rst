
IDmut Report
------------

IDmut Report: Fields
^^^^^^^^^^^^^^^^^^^^

======================================================================== ====================
Name                                                                     Data Type           
======================================================================== ====================
Sample                                                                   str                 
Branches                                                                 dict[str, str]      
Reference                                                                str                 
Discard reads with mapping qualities below this threshold                int                 
Specify the Phred score encoding of FASTQ and SAM/BAM/CRAM files         int                 
Mark base calls with Phred scores lower than this threshold as ambiguous int                 
Mark each insertion on the base to its 3' (True) or 5' (False) side      bool                
Mark all ambiguous insertions and deletions (indels)                     bool                
Retain the overhangs of paired-end mates that dovetail                   bool                
Clip this many bases from the 5' end of each read                        int                 
Clip this many bases from the 3' end of each read                        int                 
Discard alignment maps with fewer than this many reads                   int                 
Number of reads in SAM/BAM/CRAM file                                     int                 
Number of reads processed successfully                                   int                 
Number of batches                                                        int                 
Checksums of batches (SHA-512)                                           dict[str, list[str]]
Checksum of reference sequence (SHA-512)                                 str                 
Time began                                                               str                 
Time ended                                                               str                 
Time taken (minutes)                                                     float               
Version of SEISMIC-RNA                                                   str                 
======================================================================== ====================

IDmut Report: Example
^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": "",
            "idmut": ""
        },
        "Reference": "myref",
        "Discard reads with mapping qualities below this threshold": 25,
        "Specify the Phred score encoding of FASTQ and SAM/BAM/CRAM files": 33,
        "Mark base calls with Phred scores lower than this threshold as ambiguous": 25,
        "Mark each insertion on the base to its 3' (True) or 5' (False) side": true,
        "Mark all ambiguous insertions and deletions (indels)": true,
        "Retain the overhangs of paired-end mates that dovetail": true,
        "Clip this many bases from the 5' end of each read": 4,
        "Clip this many bases from the 3' end of each read": 4,
        "Discard alignment maps with fewer than this many reads": 1000,
        "Number of reads in SAM/BAM/CRAM file": 3595,
        "Number of reads processed successfully": 3595,
        "Number of batches": 1,
        "Checksums of batches (SHA-512)": {
            "idmut": [
                "8726f5b45a8588b5d462e47be651f232343214f7f755d816f6a2803e3b6bbdf17bb808ab86a74444ca2c6fcb7ef8a97f14482824b6aa524747f23fe73ab18635"
            ]
        },
        "Checksum of reference sequence (SHA-512)": "0686c09f78bae859b81f5527501018fe9768ef7595cbe1ed966aa549e227daeef8a85bf3a3dd01b9d487286fcf319e0c71ee25516d3d16b10acfe72e6eb7bcd7",
        "Time began": "2026-05-31 at 12:53:11",
        "Time ended": "2026-05-31 at 12:53:12",
        "Time taken (minutes)": 0.01,
        "Version of SEISMIC-RNA": "0.25.3"
    }
