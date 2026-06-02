********************************************************************************
Why SEISMIC-RNA
********************************************************************************

We built SEISMIC-RNA because none of the existing software for mutational
profiling data met all our needs.
A piece of software compatible with Linux or macOS that could analyze multiple
samples (FASTQ files) from many transcripts simultaneously, in parallel.
One that would perform all steps of data processing, from alignment to structure
prediction, autonomously with a single command.
And that would correct known biases in DMS-MaPseq and SHAPE-MaP data, applying
rigorous filters to ensure only high-quality input data wound up in the results.

From these requirements, we designed SEISMIC-RNA to be both easy to use yet
highly customizable, and both speedy yet accurate.
You can install it via pip and conda, then run it with a single command
(``seismic wf``) that performs all steps from alignment to structure prediction
using default settings that work in most cases.
For more advanced usage, you can run steps individually and customize dozens of
options to best suit your data.
And if the command-line interface doesn't do everything you need, SEISMIC-RNA
also provides a Python API (``import seismicrna``) for access to all its
classes and functions.
We hope that by having a low floor and high ceiling, SEISMIC-RNA will make it
easy to investigate the structure ensembles into which your favorite RNAs fold.


What SEISMIC-RNA does
================================================================================

Given sequencing reads from a chemical probing experiment (such as DMS-MaPseq
or SHAPE-MaP) and one or more reference sequences, SEISMIC-RNA:

0. Optionally demultiplexes the FASTQ files.
1. Aligns the reads to the references.
2. Identifies the mutations in each read.
3. Filters reads and positions on user-specified criteria.
4. Clusters reads by their mutations to infer structure ensembles.
5. Predicts secondary structure(s) using the mutation rates.
6. Generates a variety of graphs and RNA structure diagrams.

The same pipeline runs from a single :doc:`wf </use/workflow/wf>` command or
step-by-step from any intermediate point.


What distinguishes SEISMIC-RNA
================================================================================

End-to-end and reproducible
    The entire workflow above can be run with one command, ``seismic wf``,
    removing the need to write a script with multiple commands.
    Every step saves its output files so that the workflow can be resumed from
    the middle if needed.
  
Scales to many samples and transcripts
    SEISMIC-RNA searches directories recursively for input files, allowing you
    to provide many input files at once.
    It then processes input files in parallel to maximize speed.
  
Automated organization
    SEISMIC-RNA automatically organizes files using the names of the given
    samples and transcripts, sidestepping the need for manual organization.

Structure ensemble inference
    SEISMIC-RNA's clustering step determines whether the RNA folds into multiple
    structures, which is essential for modeling RNA structures accurately.

Scales to long RNAs
    SEISMIC-RNA can handle both short amplicons and transcripts tens of
    thousands of nucleotides long.

Open-source and scriptable
    SEISMIC-RNA is released under the GPL (see :doc:`/about/license`) and
    exposes both a CLI (see :doc:`/cli`) and a Python API (see
    :doc:`/api/index`).


Where to go next
================================================================================

- New to SEISMIC-RNA? Start with :doc:`/install/index`, then work through a
  tutorial in :doc:`/tutorials/index`.
- Looking for a specific command's options or error messages?
  See :doc:`/use/index`.
- Wondering how the science works? See :doc:`/works/index` for a visual tour
  and :doc:`/algos/index` for the algorithmic details.
- Citing SEISMIC-RNA in a paper? See :doc:`/about/cite`.
