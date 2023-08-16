
Examples
++++++++

CLI
---

Basic usage showing how to create mutation vectors from the BAM file ``DMS_100mM/MAPT-3R.bam``, with the reference sequence defined in ``human_transcriptome.fasta``::

    dreem vector --bamf DMS_100mM/MAPT-3R.bam --fasta human_transcriptome.fasta

Basic usage showing how to create mutation vectors from every BAM file in the directory ``DMS_200mM``::

    dreem vector --bamd DMS_200mM --fasta human_transcriptome.fasta

Python
------

Basic usage of the Python API to perform the same operations as shown in the CLI example. Note that ``bamf`` and ``bamd`` must each be a ``tuple`` of ``str``, not a ``str``.

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import relate
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mutvec
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import mut
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")

>>> from dreem import vector
>>> report_files_100 = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                    fasta="human_transcriptome.fasta")
>>> report_files_200 = vector.main.run(bamd=("DMS_200mM",),
...                                    fasta="human_transcriptome.fasta")
