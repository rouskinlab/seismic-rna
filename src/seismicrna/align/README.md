# DREEM Alignment Module
Contributor: Matty Allan, Yves Martin

## Purpose
Convert sequencing files into BAM files. 

Each read is compared and aligned to the reference sequence of the fasta file. Results are stored in one BAM file per reference. Additional information such as the fastqc report and the .bam/bai file is also outputed. 

## Interface

### Input Files
- [=1] ```my_fasta.fasta```. Sequence record file that contains the names and sequences of all references.
- [=1] ```my_fastq_R1.fastq```. Sequence alignment file(s) containing one or several sequences. 
- [=1] ```my_fastq_R2.fastq```. Sequence alignment file(s) containing one or several sequences. 

### Output Files

[=1] `/{sample}` Sequence alignment map file(s) folder. `{reference_k}` are the references of the fasta file found in the fastq files. 

```bash
{out_dir}:= path/to/{sample}/
  |- reference_1.bam
  |- reference_1.bam.bai
  |- reference_1_fastqc_report.txt
  |- reference_2.bam
  |- reference_2.bam.bai
  |- reference_2_fastqc_report.txt
  |- ...
```

### Command-line usage

```dreem-alignment --fastq1 [path to file] --fastq2 [path to file] --fasta [path to file]  --out_dir [path]```

- ```dreem-alignment```: Wrapper for ```run``` function in ```dreem/alignment/run.py```. 
- [=1] `--fasta / -fa` : ```my_fasta.fasta```
- [=1] `--fastq1 / -fq1`: ```my_fastq_R1.fastq```
- [=1] `--fastq2 / -fq2`: ```my_fastq_R2.fastq```
- [≤1] `--out_dir / -o`: output repository. Last directory in the path is the sample name.
- [≤1] `--fastqc_thresh`: filter out reads whose fastqc scorse do not reach that threashold.
- `[bowtie2 args]`
- `[cutadapt args]`
