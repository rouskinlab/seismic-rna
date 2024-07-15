
Align Report
------------------------------------------------------------------------

Align Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

==================================================================== ==============
Name                                                                 Data Type
==================================================================== ==============
Name of Sample                                                       str
Use demultiplexed mode                                               bool
Use paired-end mode                                                  bool
Phred score encoding in FASTQ and SAM/BAM/CRAM files                 int
Run FastQC on the initial and trimmed FASTQ files                    bool
Trim reads with Cutadapt before alignment                            bool
Phred score for read 1 quality trimming with Cutadapt                int
Phred score for read 2 quality trimming with Cutadapt                int
5' adapter for read 1 adapter trimming with Cutadapt                 list[str]
3' adapter for read 1 adapter trimming with Cutadapt                 list[str]
5' adapter for read 2 adapter trimming with Cutadapt                 list[str]
3' adapter for read 2 adapter trimming with Cutadapt                 list[str]
Minimum overlap of read and adapter during trimming with Cutadapt    int
Error tolerance for adapters during trimming with Cutadapt           float
Allow indels in adapters during trimming with Cutadapt               bool
Trim high-quality Gs from 3' end during trimming with Cutadapt       bool
Discard reads in which an adapter was found by Cutadapt              bool
Discard reads in which no adapter was found by Cutadapt              bool
Minimum length of a read to keep it after trimming with Cutadapt     int
Run Bowtie2 in local mode                                            bool
Output discordant alignments from Bowtie2                            bool
Attempt to align individual mates of unaligned pairs with Bowtie2    bool
Treat dovetailed mate pairs as concordant with Bowtie2               bool
Treat nested mate pairs as concordant with Bowtie2                   bool
Minimum score for a valid alignment with Bowtie2                     str
Minimum fragment length for valid paired-end alignments with Bowtie2 int
Maximum fragment length for valid paired-end alignments with Bowtie2 int
Minimum distance of an indel from the end of a read with Bowtie2     int
Seed length for Bowtie2                                              int
Seed interval for Bowtie2                                            str
Maximum number of consecutive failed seed extensions with Bowtie2    int
Maximum number of re-seeding attempts with Bowtie2                   int
Width of padding on alignment matrix (to allow indels) with Bowtie2  int
Valid orientations of paired-end mates with Bowtie2                  str
Output unaligned reads from Bowtie2 to a FASTQ file                  bool
Minimum mapping quality to keep an aligned read from Bowtie2         int
Number of reads initially                                            int
Number of reads after trimming                                       int
Number of reads after alignment                                      dict[str, int]
Number of reads after filtering                                      dict[str, int]
Number of reads aligned by reference                                 dict[str, int]
Branches                                                             list[str]
Time Began                                                           str
Time Ended                                                           str
Time Taken (minutes)                                                 float
Version of SEISMIC-RNA                                               str
==================================================================== ==============

Align Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Name of Sample": "ldi",
        "Use demultiplexed mode": false,
        "Use paired-end mode": true,
        "Phred score encoding in FASTQ and SAM/BAM/CRAM files": 33,
        "Run FastQC on the initial and trimmed FASTQ files": true,
        "Trim reads with Cutadapt before alignment": true,
        "Phred score for read 1 quality trimming with Cutadapt": 25,
        "Phred score for read 2 quality trimming with Cutadapt": 25,
        "5' adapter for read 1 adapter trimming with Cutadapt": [
            "GCTCTTCCGATCT"
        ],
        "3' adapter for read 1 adapter trimming with Cutadapt": [
            "AGATCGGAAGAGC"
        ],
        "5' adapter for read 2 adapter trimming with Cutadapt": [
            "GCTCTTCCGATCT"
        ],
        "3' adapter for read 2 adapter trimming with Cutadapt": [
            "AGATCGGAAGAGC"
        ],
        "Minimum overlap of read and adapter during trimming with Cutadapt": 6,
        "Error tolerance for adapters during trimming with Cutadapt": 0.1,
        "Allow indels in adapters during trimming with Cutadapt": true,
        "Trim high-quality Gs from 3' end during trimming with Cutadapt": false,
        "Discard reads in which an adapter was found by Cutadapt": false,
        "Discard reads in which no adapter was found by Cutadapt": false,
        "Minimum length of a read to keep it after trimming with Cutadapt": 20,
        "Run Bowtie2 in local mode": true,
        "Output discordant alignments from Bowtie2": false,
        "Attempt to align individual mates of unaligned pairs with Bowtie2": false,
        "Treat dovetailed mate pairs as concordant with Bowtie2": false,
        "Treat nested mate pairs as concordant with Bowtie2": true,
        "Minimum score for a valid alignment with Bowtie2": "L,1,0.5",
        "Minimum fragment length for valid paired-end alignments with Bowtie2": 0,
        "Maximum fragment length for valid paired-end alignments with Bowtie2": 600,
        "Minimum distance of an indel from the end of a read with Bowtie2": 4,
        "Seed length for Bowtie2": 20,
        "Seed interval for Bowtie2": "L,1,0.1",
        "Maximum number of consecutive failed seed extensions with Bowtie2": 4,
        "Maximum number of re-seeding attempts with Bowtie2": 2,
        "Width of padding on alignment matrix (to allow indels) with Bowtie2": 2,
        "Valid orientations of paired-end mates with Bowtie2": "fr",
        "Output unaligned reads from Bowtie2 to a FASTQ file": true,
        "Minimum mapping quality to keep an aligned read from Bowtie2": 25,
        "Number of reads initially": 566520,
        "Number of reads after trimming": 565068,
        "Number of reads after alignment": {
            "reads, were paired": 565068,
            "reads, were paired, aligned concordantly 0 times": 2613,
            "reads, were paired, aligned concordantly exactly 1 time": 562448,
            "reads, were paired, aligned concordantly >1 times": 7
        },
        "Number of reads after filtering": {
            "paired-end, both mates mapped": 562325,
            "paired-end, one mate unmapped": 0
        },
        "Number of reads aligned by reference": {
            "sars2_1799": 562325
        },
        "Branches": [],
        "Time Began": "2023-11-10 at 20:16:29",
        "Time Ended": "2023-11-10 at 20:17:39",
        "Time Taken (minutes)": 1.16,
        "Version of SEISMIC-RNA": "0.9.3"
    }
