

Vectoring Report
++++++++++++++++

The vectoring report file summarizes the results of vectoring and enables DREEM to locate all mutation vector batch files and verify their integrity.
The report contains the following fields:

- ``Sample``: name of the sample from which the reads came
- ``Reference``: name of the reference sequence to which the reads were aligned
- ``Length (nt)``: number of nucleotides in the reference sequence
- ``Sequence``: full sequence of the reference
- ``Vectors``: number of reads (or pairs of paired-end reads) that were successfully converted into mutation vectors and written to an ORC file
- ``Errors``: number of reads (or pairs of paired-end reads) that failed to be converted into mutation vectors
- ``Percent``: percentage of total reads that were vectorized successfully: ``"Reads Vectorized" / ("Reads Vectorized" + "Reads with Errors")``
- ``Batches``: number of ORC files written, each containing a batch of mutation vectors 
- ``MD5 Checksums``: comma-separated list of the MD5 checksum of each ORC file (for verifying file integrity)
- ``Began``: local time at which vectoring began, formatted as ``YYYY-MM-DD hh:mm:ss.µs``
- ``Ended``: local time at which vectoring ended, formatted as ``YYYY-MM-DD hh:mm:ss.µs``
- ``Time (s)``: time taken to complete vectorization, in seconds: ``"Ended" - "Began"``
- ``Speed (1/s)``: average number of vectors written per second: ``"Reads Vectorized" / "Time Taken (s)"``

**Example: report.json**

Example contents of a vectoring report file::

    {
     "Sample": "DMS_100mM",
     "Reference": "SARS2",
     "Length (nt)": "30",
     "Sequence": "AGTCTGTACCGTCTGCGGTATGTGGAAAGG",
     "Vectorized": "272638",
     "Errors": "0",
     "Percent": "100.0",
     "Batches": "2",
     "MD5 Checksums": "9f4ffa03d172df618fa664358cc9dc8e, a8f23d7cf38e8a15bcb298037bd23f5c",
     "Began": "2023-03-22 20:00:05.712658",
     "Ended": "2023-03-22 20:00:28.427814",
     "Time (s)": "22.715",
     "Speed (1/s)": "12002.6"
    }
