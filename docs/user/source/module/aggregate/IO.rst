
I/O files
++++++++++++++++++++++++

*This module needs the files to be sorted into a specific folder structure. The folder structure is described below.*

**Input files structure**

As input, aggregate requires one or more files of :ref:`bitvector` batches for each section of a reference. The folder structure is as follows::

    sample_1/              # <=> a fastq file <<< give this path to aggregate the entire sample
        |- reference_1/    # <=> a single reference in the fasta file
            |- section_1/  # <=> a section (a sub-sequence of a reference) 
                |- 0.orc   # <=> a batch of bitvectors
                |- 1.orc
                |- report.json  # <=> a report file <<< give this path to aggregate a single section
            |- section_2  # <=> another section for this reference
                |- 0.orc
                |- report.json
        |- reference_2    # <=> another reference
            |- section_3
                |- 0.orc
                |- 1.orc
                |- 2.orc
                |- report.json


**Output file structure**

:: 

    sample_1.json
