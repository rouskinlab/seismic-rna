
Examples
++++++++

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