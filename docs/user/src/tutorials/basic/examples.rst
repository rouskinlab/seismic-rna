Examples
++++++++

CLI
---

Basic usage of the align module, showing how to align a pair of FASTQ files containing paired-end reads
(``DMS_100mM_R1.fq.gz`` and ``DMS_100mM_R2.fq.gz``) to a set of reference sequences in ``human_transcriptome.fasta``::

    dreem align --fastq1 DMS_100mM_R1.fq.gz --fastq2 DMS_100mM_R2.fq.gz --fasta human_transcriptome.fasta

More advanced usage, showing how to align every pair of FASTQ files in a directory of demultiplexed mate 1 and mate 2 FASTQ files
(``DMS_200mM``) to a set of reference sequences in ``human_transcriptome.fasta``::

    dreem align --fastq12_dir DMS_200mM --fasta human_transcriptome.fasta


Python
------

Basic usage of the Python API, showing the same examples as with the CLI:

>>> from dreem import align
>>> bam_files_100 = align.main.run(fastqm=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(dmfastqm=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")m=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(dmfastqm=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")

>>> from dreem import align
>>> bam_files_100 = align.main.run(fastq1=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(dmfastqm=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")

>>> from dreem import align
>>> bam_files_100 = align.main.run(fastq1=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(fastq12_dir=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")

>>> from dreem import align
>>> bam_files_100 = align.main.run(fastq1=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(fastq12_dir=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")

>>> from dreem import align
>>> bam_files_100 = align.main.run(fastqm=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(fastq12_dir=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")

>>> from dreem import align
>>> bam_files_100 = align.main.run(fastqm=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(fastq12_dir=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")

>>> from dreem import align
>>> bam_files_100 = align.main.run(fastq1=("DMS_100mM_R1.fq.gz",),
...                                fastq2=("DMS_100mM_R2.fq.gz",),
...                                fasta="human_transcriptome.fasta")
>>> bam_files_200 = align.main.run(fastq12_dir=("DMS_200mM",),
...                                fasta="human_transcriptome.fasta")
