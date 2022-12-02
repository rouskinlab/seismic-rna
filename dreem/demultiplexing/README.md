# DREEM Demultiplexing Module
Contributor: Scott Grote, Yves Martin

## Purpose
- Group the reads by barcode

## Interface

### Input Files
- [≥1] ```{sample_x}_R1.fastq```. Forward primer sequence alignment file(s) generated by the sequencer. There must be one file for each sample, and its name must be ```{sample}_R1.fastq```.  
- [≥1] ```{sample_x}_R2.fastq```. Forward primer sequence alignment file(s) generated by the sequencer. There must be one file for each sample, and its name must be ```{sample}_R2.fastq```.  
- [=1] ```library.csv```. A CSV file with the following columns: 
   - `construct`: name of the output fasta file.
   - `barcode_start`: a 0-index of the beginning of the barcode sequence.
   - `barcode_end`: a 0-index of the end of the barcode sequence (this index being non-included).
   - `barcode`: a string of A C G T forming the barcode.

### Output files
- [≥1] `{sample_x}/{construct_xy}_R1.fastq`. Sequence alignment file(s) containing reads from `{sample_x}_R1.fastq` with `barcode` of row `construct_xy` as a barcode.
- [≥1] `{sample_x}/{construct_xy}_R2.fastq`. Sequence alignment file(s) containing reads from `{sample_x}_R2.fastq` with `barcode` of row `construct_xy` as a barcode.

### Command-line usage

```dreem-demultiplexing —-fastq1 [file] --fastq2 [file] —-library [file] --out_dir [sample]```

- ```dreem-demultiplexing```: Wrapper for ```run``` function in ```dreem/demultiplexing/run.py```. 
- [≥1] `--fastq1`: ```{sample_x}_R1.fastq```
- [≥1] `--fastq2`: ```{sample_x}_R2.fastq```
- [=1] `--library` : ```library.csv```
- [≤1] `--out_dir`: Name of the output directory.
- [≤1] `--barcode_start`: Start position of the barcode in the read (uncompatible with library)
- [≤1] `--barcode_end`: End position of the barcode in the read (uncompatible with library)