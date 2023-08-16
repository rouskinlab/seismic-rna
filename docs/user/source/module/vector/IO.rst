
I/O files
+++++++++

*This module needs the files to be sorted into a specific folder structure. The folder structure is described below.*

**Input file structure**

The module takes one or more files of aligned reads in binary alignment map (BAM) format.
Each file must be named after the reference sequence to which it was aligned, end with the file extension ``.bam``, and reside in a directory named after the sample from which the reads came::

    |- sample_1/
        |- reference_1.bam
        |- reference_2.bam
    |- sample_2/
        |- reference_1.bam
        |- reference_3.bam


**Output file structure**

For each sample aligned to each reference, zero or more mutation vector batch files in `Apache ORC format <https://orc.apache.org/docs>`_ (numbered ``1.orc``, ``2.orc``, ...) are written.
One report file (``report.json``) in `JSON format <https://www.json.org/json-en.html>`_) is also written and contains metadata about the vectors.
The following example illustrates a possible output file structure::

    |- sample_1/
        |- reference_1/
            |- report.json
            |- 1.orc
            |- 2.orc
        |- reference_2
            |- report.json
            |- 1.orc
            |- 2.orc
            |- 3.orc
    |- sample_2/
        |- reference_2/
            |- report.json
            |- 1.orc
        |- reference_3
            |- report.json
            |- 1.orc
            |- 2.orc
