
I/O files
++++++++++++++++++++++++

*This module needs the files to be sorted into a specific folder structure. The folder structure is described below.*

**Input**


As input, demultiplex requires either a library.csv file with the barcode_start_index and barcode_length, but this can be alternatively be given as input in the CLI::

    sample_1_R1              # <=> a fastq file <<< give this path to be demultiplexed
    sample_1_R2              # <=> a fastq file <<< give this path to be demultiplexed
    sample_ref               # <=> a fasta file <<< give this path of references 
        

**Output**

As output, demultiplex provides a lot of temporary files but the primary output is the demultiplexed fastqs. The folder structure is as follows::

    temp/              # <=> a fastq file <<< give this path to aggregate the entire sample
        |- sample_1/    # <=> holds all the demultiplex temporary files
            |- sample_1_R1_pickle_split/  # <=> holds paritioned fastqs stored as dictionaries, usually split into 10
                |- split_0.fastq   # <=> a split fastq
                |- split_0.p       # <=> a split fastq stored as dictionary as python pickle
                |- split_1.fastq   
                |- split_1.p
            |- sample_1_R2_pickle_split/  
                |- split_0.fastq   
                |- split_0.p       
                |- split_1.fastq   
                |- split_1.p
            |- sample_fqs/    # <=> holds all 
                |- reference_1_R1.fastq      # <=> fastq representing reads found from sample_1_R1
                |- reference_1_R2.fastq      # <=> fastq representing reads found from sample_1_R2
                |- reference_2_R1.fastq
                |- reference_2_R2.fastq
            |- sequence_data/    # <=> holds all the temporary data used in demultiplexing
                |- reference_1/     # <=> the references temporary data
                    |- fq1/
                        |- specific_pattern_grepped.fastq # <=> demultiplexing requires several different types of grepping to isolate what reads belong to what reference
                        |- complete_set_of_reads.p   # <=> the total set of read ids found post 
                        |- read_id_data.p  # <=> the read ids found by what type of grep
                        |- unfiltered.fastq # <=> the complete set of reads before a secondary signiture filtering process 
                    |- fq2/
                        |- specific_pattern_grepped.fastq
                        |- complete_set_of_reads.p
                        |- read_id_data.p 
                        |- unfiltered.fastq
                |- reference_2/
                    |- fq1/
                        |- specific_pattern_grepped.fastq
                        |- complete_set_of_reads.p
                        |- read_id_data.p 
                        |- unfiltered.fastq
                    |- fq2/
                        |- specific_pattern_grepped.fastq
                        |- complete_set_of_reads.p
                        |- read_id_data.p 
                        |- unfiltered.fastq
                |- multigrepped_reads.p     # <=> dictionary mapping read ids to a list of the references to which it mapped
            |-demultiplex_info.csv          # <=> a csv showing the read counts found for each type of query

