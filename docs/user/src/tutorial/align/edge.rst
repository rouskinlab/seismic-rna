
Edge cases handling
++++++++++++++++++++++++

Zero FASTQ files
----------------
The alignment module will complete without crashing (albeit with several warning/error messages) if no FASTQ files are given.

Empty FASTQ File
----------------
DREEM can handle FASTQ files with no reads because every piece of software that the alignment module uses on the backend
(namely FASTQC, Cutadapt, Bowtie2, and Samtools) can handle empty FASTQ files.
One algorithm (for deduplication) is implemented in DREEM itself; this algorithm can also handle empty files, but will crash if given a corrupted or misformatted file.
All expected output files will be generated, but they will contain no data (e.g. only SAM header information).
