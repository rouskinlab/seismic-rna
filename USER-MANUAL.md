# SEISMIC-RNA User Manual


## List of steps in the pipeline

![Flowchart of SEISMIC-RNA](flowchart.png "Flowchart")

1. `align`: (optional) Align sequencing reads to one or more references, after optional quality control and trimming.
2. `relate`: Determine the relationship (e.g. match, substitution, not covered) between each position in the reference and the base call to which it aligned in each read.
3. `mask`: Specify which relationships count as mutations, then remove reads and positions that fail filters based on mutations and coverage.
4. `cluster`: (optional) Cluster the mutations to identify alternative RNA structures.
5. `table`: Output two tables counting each type of relationship 1) per position and 2) per read.
6. `struct`: (optional) Predict the secondary structure of the RNA given the mutation rates.
7. `graph`: (optional) Graph data from tables.


## General tips

- For each step, you can display a list of the required arguments and all available options by typing the command followed by `--help`: for example, `seismic all --help`.


## Run all steps automatically


### Command
```commandline
seismic all refs.fa -x paired_R1.fq.gz -x paired_R2.fq.gz
```
or
```commandline
seismic all refs.fa -y interleaved.fq.gz
```
or
```commandline
seismic all refs.fa -z single.fq.gz
```

Any option for one of the individual steps (described below) can be given as an option to `all`.


## Run the step `align`

### Command
```commandline
seismic all refs.fa -x paired_R1.fq.gz -x paired_R2.fq.gz
```
or
```commandline
seismic all refs.fa -y interleaved.fq.gz
```
or
```commandline
seismic all refs.fa -z single.fq.gz
```

### Descriptions
- `refs.fa` is a file of all reference sequences in FASTA format.
- `paired_R1.fq.gz` is a file of paired-end first reads in (possibly gzipped) FASTQ format.
- `paired_R2.fq.gz` is a file of paired-end second reads in (possibly gzipped) FASTQ format.
- `interleaved.fq.gz` is a file of interleaved paired-end reads in (possibly gzipped) FASTQ format.
- `single.fq.gz` is a file of single-end reads in (possibly gzipped) FASTQ format.


### Tips

- You can align multiple FASTQ files (or pairs of paired-end FASTQ files) simultaneously by
  - giving a directory, which will be searched recursively for FASTQ files
  - giving `-x`, `-y`, and/or `-z` multiple times, each with a different file/directory
  - a combination of the above
- For a FASTQ file given via `-y` or `-z`, the name of the file up to the file extension becomes the sample name.
- For a FASTQ file given via `-x`, the mate number (1 or 2) must be specified as `_R1`/`_R2`, `_mate1`/`_mate2`, or `_1_sample`/`_2_sample` immediately before the file extension. The name of the file, minus the extension and the mate specifier, becomes the sample name.
- FASTQ files tend to be very large (100s of megabytes to 100s of gigabytes). Compressing FASTQ files with `gzip` can reduce the file size to less than one third of the original and is _highly recommended_. Gzipped FASTQ files (`.fastq.gz`/`.fq.gz`) can be fed directly into SEISMIC-RNA without being decompressed.


## Run the step `relate`

### Command
```commandline
seismic relate refs.fa sample1.bam ...
```

### Descriptions
- `refs.fa` is a file of all reference sequences in FASTA format.
- `sample1.bam` is an alignment map file in BAM format.


### Tips

- You can relate multiple BAM files simultaneously by
  - giving a directory, which will be searched recursively for BAM files
  - giving multiple arguments, separated by spaces
  - using the shell's built-in globbing capabilities
  - a combination of the above


## Run the step `mask`

### Command
```commandline
seismic mask report-relate.json ...
```

### Descriptions
- `report-relate.json` is a report file output by the `relate` step.


### Tips

- You can mask multiple `relate` reports simultaneously by
  - giving a directory, which will be searched recursively for `relate` reports
  - giving multiple arguments, separated by spaces
  - using the shell's built-in globbing capabilities
  - a combination of the above


## Run the step `cluster`

### Command
```commandline
seismic cluster report-mask.json ... -k 2
```

### Descriptions
- `report-mask.json` is a report file output by the `mask` step.
- `-k` is the maximum number of clusters to attempt (default: 2).


### Tips

- You can cluster multiple `mask` reports simultaneously by
  - giving a directory, which will be searched recursively for `mask` reports
  - giving multiple arguments, separated by spaces
  - using the shell's built-in globbing capabilities
  - a combination of the above


## Run the step `table`

### Command
```commandline
seismic cluster report-*.json ...
```

### Descriptions
- `report-*.json` is a report file output by the `relate`, `mask`, or `cluster` step.


### Tips

- You can tabulate multiple reports simultaneously by
  - giving a directory, which will be searched recursively for reports
  - giving multiple arguments, separated by spaces
  - using the shell's built-in globbing capabilities
  - a combination of the above
- To reduce the running time, you can specify that only certain relationships be tabulated using the option `-r` followed by a string of one-letter codes for the relationships to tabulate: for example, `-r vrms`. The one-letter codes are as follows:
  - `v`: Covered by the read
  - `r`: Matched the reference
  - `m`: Mutated (any type)
  - `s`: Substituted (any type)
  - `a`: Substituted to A
  - `c`: Substituted to C
  - `g`: Substituted to G
  - `t`: Substituted to T
  - `d`: Deleted
  - `i`: Inserted


## Run the step `struct`

### Command
```commandline
seismic struct table.csv ... 
```

### Descriptions
- `table.csv` is a table in CSV format output by the `table` step.


### Tips

- You can predict structures from multiple tables simultaneously by
  - giving a directory, which will be searched recursively for tables
  - giving multiple tables, separated by spaces
  - using the shell's built-in globbing capabilities
  - a combination of the above


## Run the step `graph`

### Command
```commandline
seismic graph seq table.csv ... 
```

### Descriptions
- `table.csv` is a table in CSV format output by the `table` step.


### Tips

- You can graph multiple tables simultaneously by
  - giving a directory, which will be searched recursively for tables
  - giving multiple tables, separated by spaces
  - using the shell's built-in globbing capabilities
  - a combination of the above
- You can graph a subset of the relationships in the table using the option `-r` followed by a string of one-letter codes for the relationships to tabulate: for example, `-r acgtdi`. The one-letter codes are as follows:
  - `v`: Covered by the read
  - `n`: Informed (i.e. unambiguously matched or mutated)
  - `r`: Matched the reference
  - `m`: Mutated (any type)
  - `s`: Substituted (any type)
  - `a`: Substituted to A
  - `c`: Substituted to C
  - `g`: Substituted to G
  - `t`: Substituted to T
  - `d`: Deleted
  - `i`: Inserted
