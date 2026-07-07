
Align Report
------------

Align Sample Report
^^^^^^^^^^^^^^^^^^^

Align Sample Report: Fields
===========================

======================================================================================= ==============
Name                                                                                    Data Type     
======================================================================================= ==============
Sample                                                                                  str           
Branches                                                                                dict[str, str]
Seed for the random number generator                                                    int           
Use demultiplexed mode                                                                  bool          
Use paired-end mode                                                                     bool          
Specify the Phred score encoding of FASTQ and SAM/BAM/CRAM files                        int           
Use fastp to QC, filter, and trim reads before alignment                                bool          
Trim low-quality bases from the 5' ends of reads                                        bool          
Trim low-quality bases from the 3' ends of reads                                        bool          
Use this window size (nt) for --fastp-5 and --fastp-3                                   int           
Use this mean quality threshold for --fastp-5 and --fastp-3                             int           
Trim poly(G) tails (two-color sequencing artifacts) from the 3' end                     str           
Minimum number of Gs to consider a poly(G) tail for --fastp-poly-g                      int           
Trim poly(X) tails (i.e. of any nucleotide) from the 3' end                             bool          
Minimum number of bases to consider a poly(X) tail for --fastp-poly-x                   int           
Trim adapter sequences from the 3' ends of reads                                        bool          
Trim this adapter sequence from the 3' ends of read 1s                                  str           
Trim this adapter sequence from the 3' ends of read 2s                                  str           
Trim adapter sequences in this FASTA file from the 3' ends of reads                     str           
Automatically detect the adapter sequences for paired-end reads                         bool          
Discard reads shorter than this length                                                  int           
Align reads in local mode rather than end-to-end mode                                   bool          
Output paired-end reads whose mates align discordantly                                  bool          
Attempt to align individual mates of pairs that fail to align                           bool          
Consider dovetailed mate pairs to align concordantly                                    bool          
Consider nested mate pairs to align concordantly                                        bool          
Discard alignments that score below this threshold                                      str           
Discard paired-end alignments shorter than this many bases                              int           
Discard paired-end alignments longer than this many bases                               int           
Do not place gaps within this many bases from the end of a read                         int           
Use this seed length for Bowtie2                                                        int           
Seed Bowtie2 alignments at this interval                                                str           
Discard alignments if over this many consecutive seed extensions fail                   int           
Re-seed reads with repetitive seeds up to this many times                               int           
Pad the alignment matrix with this many bases (to allow gaps)                           int           
Require paired mates to have this orientation                                           str           
Output unaligned reads to a FASTQ file                                                  bool          
Discard reads with mapping qualities below this threshold                               int           
Separate each alignment map into forward- and reverse-strand reads                      bool          
With --sep-strands, consider forward mate 1s and reverse mate 2s to be forward-stranded bool          
With --sep-strands, add this label to each reverse-strand reference                     str           
Discard alignment maps with fewer than this many reads                                  int           
Number of reads in the FASTQ file(s)                                                    int           
Number of reads after trimming                                                          int           
Number of reads after alignment                                                         dict[str, int]
Number of reads after filtering                                                         dict[str, int]
Number of reads aligned to each reference                                               dict[str, int]
Checksum of the reference fasta (SHA-512)                                               str           
Checksum(s) of the input fastq(s) (SHA-512)                                             dict[str, str]
Time began                                                                              str           
Time ended                                                                              str           
Time taken (minutes)                                                                    float         
Version of SEISMIC-RNA                                                                  str           
======================================================================================= ==============

Align Sample Report: Example
============================

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": ""
        },
        "Seed for the random number generator": null,
        "Use demultiplexed mode": false,
        "Use paired-end mode": true,
        "Specify the Phred score encoding of FASTQ and SAM/BAM/CRAM files": 33,
        "Use fastp to QC, filter, and trim reads before alignment": true,
        "Trim low-quality bases from the 5' ends of reads": false,
        "Trim low-quality bases from the 3' ends of reads": true,
        "Use this window size (nt) for --fastp-5 and --fastp-3": 6,
        "Use this mean quality threshold for --fastp-5 and --fastp-3": 25,
        "Trim poly(G) tails (two-color sequencing artifacts) from the 3' end": "auto",
        "Minimum number of Gs to consider a poly(G) tail for --fastp-poly-g": 10,
        "Trim poly(X) tails (i.e. of any nucleotide) from the 3' end": false,
        "Minimum number of bases to consider a poly(X) tail for --fastp-poly-x": 10,
        "Trim adapter sequences from the 3' ends of reads": true,
        "Trim this adapter sequence from the 3' ends of read 1s": "",
        "Trim this adapter sequence from the 3' ends of read 2s": "",
        "Trim adapter sequences in this FASTA file from the 3' ends of reads": null,
        "Automatically detect the adapter sequences for paired-end reads": true,
        "Discard reads shorter than this length": 9,
        "Align reads in local mode rather than end-to-end mode": true,
        "Output paired-end reads whose mates align discordantly": false,
        "Attempt to align individual mates of pairs that fail to align": false,
        "Consider dovetailed mate pairs to align concordantly": false,
        "Consider nested mate pairs to align concordantly": true,
        "Discard alignments that score below this threshold": "L,1,0.8",
        "Discard paired-end alignments shorter than this many bases": 0,
        "Discard paired-end alignments longer than this many bases": 600,
        "Do not place gaps within this many bases from the end of a read": 4,
        "Use this seed length for Bowtie2": 20,
        "Seed Bowtie2 alignments at this interval": "L,1,0.1",
        "Discard alignments if over this many consecutive seed extensions fail": 4,
        "Re-seed reads with repetitive seeds up to this many times": 2,
        "Pad the alignment matrix with this many bases (to allow gaps)": 2,
        "Require paired mates to have this orientation": "fr",
        "Output unaligned reads to a FASTQ file": true,
        "Discard reads with mapping qualities below this threshold": 25,
        "Separate each alignment map into forward- and reverse-strand reads": false,
        "With --sep-strands, consider forward mate 1s and reverse mate 2s to be forward-stranded": false,
        "With --sep-strands, add this label to each reverse-strand reference": "-rev",
        "Discard alignment maps with fewer than this many reads": 1000,
        "Number of reads in the FASTQ file(s)": 4008,
        "Number of reads after trimming": 4008,
        "Number of reads after alignment": {
            "reads, were paired": 4008,
            "reads, were paired, aligned concordantly 0 times": 372,
            "reads, were paired, aligned concordantly exactly 1 time": 3636,
            "reads, were paired, aligned concordantly >1 times": 0
        },
        "Number of reads after filtering": {
            "paired-end, both mates mapped": 3595,
            "paired-end, one mate unmapped": 0
        },
        "Number of reads aligned to each reference": {
            "myref": 3595
        },
        "Checksum of the reference fasta (SHA-512)": "e4dcaae215ebadadc57010d22d93f4311cdc0bc02fc2f65db558a14b671d898ef8ff541817bf374aa60688627dbb9601581f7b271e60eb80cc36fc0c7dd51882",
        "Checksum(s) of the input fastq(s) (SHA-512)": {
            "fastq1": "a28f7d8cf7856603cee938c0d6a91561063d6eb8a34ebc2d9bee2b9fa8c04f9c528afe8e6b039bb62467bf71a6358d892fc5fdd2bec9e51902a116f67b791a9d",
            "fastq2": "15242a00abef156c10b18bef6717cd2e9b1ee54165d89aa490ea0176e7977d7f24a3c9e12dece493c67a6ab77b9a43f7ad291a60189275266d22a31882ef1046"
        },
        "Time began": "2026-05-31 at 12:53:08",
        "Time ended": "2026-05-31 at 12:53:09",
        "Time taken (minutes)": 0.01,
        "Version of SEISMIC-RNA": "0.25.3"
    }

Align Reference Report
^^^^^^^^^^^^^^^^^^^^^^

Align Reference Report: Fields
==============================

======================================================================================= ==============
Name                                                                                    Data Type     
======================================================================================= ==============
Sample                                                                                  str           
Branches                                                                                dict[str, str]
Reference                                                                               str           
Seed for the random number generator                                                    int           
Use demultiplexed mode                                                                  bool          
Use paired-end mode                                                                     bool          
Specify the Phred score encoding of FASTQ and SAM/BAM/CRAM files                        int           
Use fastp to QC, filter, and trim reads before alignment                                bool          
Trim low-quality bases from the 5' ends of reads                                        bool          
Trim low-quality bases from the 3' ends of reads                                        bool          
Use this window size (nt) for --fastp-5 and --fastp-3                                   int           
Use this mean quality threshold for --fastp-5 and --fastp-3                             int           
Trim poly(G) tails (two-color sequencing artifacts) from the 3' end                     str           
Minimum number of Gs to consider a poly(G) tail for --fastp-poly-g                      int           
Trim poly(X) tails (i.e. of any nucleotide) from the 3' end                             bool          
Minimum number of bases to consider a poly(X) tail for --fastp-poly-x                   int           
Trim adapter sequences from the 3' ends of reads                                        bool          
Trim this adapter sequence from the 3' ends of read 1s                                  str           
Trim this adapter sequence from the 3' ends of read 2s                                  str           
Trim adapter sequences in this FASTA file from the 3' ends of reads                     str           
Automatically detect the adapter sequences for paired-end reads                         bool          
Discard reads shorter than this length                                                  int           
Align reads in local mode rather than end-to-end mode                                   bool          
Output paired-end reads whose mates align discordantly                                  bool          
Attempt to align individual mates of pairs that fail to align                           bool          
Consider dovetailed mate pairs to align concordantly                                    bool          
Consider nested mate pairs to align concordantly                                        bool          
Discard alignments that score below this threshold                                      str           
Discard paired-end alignments shorter than this many bases                              int           
Discard paired-end alignments longer than this many bases                               int           
Do not place gaps within this many bases from the end of a read                         int           
Use this seed length for Bowtie2                                                        int           
Seed Bowtie2 alignments at this interval                                                str           
Discard alignments if over this many consecutive seed extensions fail                   int           
Re-seed reads with repetitive seeds up to this many times                               int           
Pad the alignment matrix with this many bases (to allow gaps)                           int           
Require paired mates to have this orientation                                           str           
Output unaligned reads to a FASTQ file                                                  bool          
Discard reads with mapping qualities below this threshold                               int           
Separate each alignment map into forward- and reverse-strand reads                      bool          
With --sep-strands, consider forward mate 1s and reverse mate 2s to be forward-stranded bool          
With --sep-strands, add this label to each reverse-strand reference                     str           
Discard alignment maps with fewer than this many reads                                  int           
Number of reads in the FASTQ file(s)                                                    int           
Number of reads after trimming                                                          int           
Number of reads after alignment                                                         dict[str, int]
Number of reads after filtering                                                         dict[str, int]
Number of reads aligned to each reference                                               dict[str, int]
Checksum of the reference fasta (SHA-512)                                               str           
Checksum(s) of the input fastq(s) (SHA-512)                                             dict[str, str]
Time began                                                                              str           
Time ended                                                                              str           
Time taken (minutes)                                                                    float         
Version of SEISMIC-RNA                                                                  str           
======================================================================================= ==============

Align Reference Report: Example
===============================

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": ""
        },
        "Reference": "myref",
        "Seed for the random number generator": null,
        "Use demultiplexed mode": true,
        "Use paired-end mode": true,
        "Specify the Phred score encoding of FASTQ and SAM/BAM/CRAM files": 33,
        "Use fastp to QC, filter, and trim reads before alignment": true,
        "Trim low-quality bases from the 5' ends of reads": false,
        "Trim low-quality bases from the 3' ends of reads": true,
        "Use this window size (nt) for --fastp-5 and --fastp-3": 6,
        "Use this mean quality threshold for --fastp-5 and --fastp-3": 25,
        "Trim poly(G) tails (two-color sequencing artifacts) from the 3' end": "auto",
        "Minimum number of Gs to consider a poly(G) tail for --fastp-poly-g": 10,
        "Trim poly(X) tails (i.e. of any nucleotide) from the 3' end": false,
        "Minimum number of bases to consider a poly(X) tail for --fastp-poly-x": 10,
        "Trim adapter sequences from the 3' ends of reads": true,
        "Trim this adapter sequence from the 3' ends of read 1s": "",
        "Trim this adapter sequence from the 3' ends of read 2s": "",
        "Trim adapter sequences in this FASTA file from the 3' ends of reads": null,
        "Automatically detect the adapter sequences for paired-end reads": true,
        "Discard reads shorter than this length": 9,
        "Align reads in local mode rather than end-to-end mode": true,
        "Output paired-end reads whose mates align discordantly": false,
        "Attempt to align individual mates of pairs that fail to align": false,
        "Consider dovetailed mate pairs to align concordantly": false,
        "Consider nested mate pairs to align concordantly": true,
        "Discard alignments that score below this threshold": "L,1,0.8",
        "Discard paired-end alignments shorter than this many bases": 0,
        "Discard paired-end alignments longer than this many bases": 600,
        "Do not place gaps within this many bases from the end of a read": 4,
        "Use this seed length for Bowtie2": 20,
        "Seed Bowtie2 alignments at this interval": "L,1,0.1",
        "Discard alignments if over this many consecutive seed extensions fail": 4,
        "Re-seed reads with repetitive seeds up to this many times": 2,
        "Pad the alignment matrix with this many bases (to allow gaps)": 2,
        "Require paired mates to have this orientation": "fr",
        "Output unaligned reads to a FASTQ file": true,
        "Discard reads with mapping qualities below this threshold": 25,
        "Separate each alignment map into forward- and reverse-strand reads": false,
        "With --sep-strands, consider forward mate 1s and reverse mate 2s to be forward-stranded": false,
        "With --sep-strands, add this label to each reverse-strand reference": "-rev",
        "Discard alignment maps with fewer than this many reads": 1000,
        "Number of reads in the FASTQ file(s)": 4008,
        "Number of reads after trimming": 4008,
        "Number of reads after alignment": {
            "reads, were paired": 4008,
            "reads, were paired, aligned concordantly 0 times": 372,
            "reads, were paired, aligned concordantly exactly 1 time": 3636,
            "reads, were paired, aligned concordantly >1 times": 0
        },
        "Number of reads after filtering": {
            "paired-end, both mates mapped": 3595,
            "paired-end, one mate unmapped": 0
        },
        "Number of reads aligned to each reference": {
            "myref": 3595
        },
        "Checksum of the reference fasta (SHA-512)": "e4dcaae215ebadadc57010d22d93f4311cdc0bc02fc2f65db558a14b671d898ef8ff541817bf374aa60688627dbb9601581f7b271e60eb80cc36fc0c7dd51882",
        "Checksum(s) of the input fastq(s) (SHA-512)": {
            "fastq1": "a28f7d8cf7856603cee938c0d6a91561063d6eb8a34ebc2d9bee2b9fa8c04f9c528afe8e6b039bb62467bf71a6358d892fc5fdd2bec9e51902a116f67b791a9d",
            "fastq2": "15242a00abef156c10b18bef6717cd2e9b1ee54165d89aa490ea0176e7977d7f24a3c9e12dece493c67a6ab77b9a43f7ad291a60189275266d22a31882ef1046"
        },
        "Time began": "2026-05-31 at 12:53:08",
        "Time ended": "2026-05-31 at 12:53:09",
        "Time taken (minutes)": 0.01,
        "Version of SEISMIC-RNA": "0.25.3"
    }
