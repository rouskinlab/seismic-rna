
Get started
++++++++++++++++


Let's aggregate ``sample_1``.

You need:
    - the paths to fastqs 
    - the reference genome under the `fasta format <https://en.wikipedia.org/wiki/FASTA_format>`_ ``reference.fasta``.
    - either: 
    - a. the library csv with a reference field, a barcode start index field, a barcode length, and optionally a secondary signiture start index and a secondary signiture length 
    - b. give barcode start index and barcode length
    - 

You get:
    -  ``sample_1.json`` under the :ref:`dreem_output` format.



CLI
---------

Run the command:

::
    
    dreem demultiplex —-fasta path/to/reference.fasta --fastq1 path/to/sample_1_R1.fastq --fastq2 path/to/sample_1_R2.fastq --barcode-start X --barcode-length Y 



Python
------------

Run the following code:

:: 

    from ..demultiplex.demultiplex import demultiplex_run
    demultiplex_run(
                    mixed_fastq1 = path/to/sample_1_R1.fastq,
                    mixed_fastq2 = path/to/sample_1_R1.fastq,
                    barcode_start = X,
                    barcode_length = Y,
                    fasta = path/to/reference.fasta)


See the :ref:`reference` for more details.
