
Input file patterns
========================================================================

For commands that take a list of input files as positional arguments,
these files can be given in a variety of ways so that you can choose a
convenient manner of specifying input files.

Example output directory
------------------------------------------------------------------------

Assume that the output directory has these contents::

    sample_1/
        align/
            ref_1.bam
            ref_2.bam
        relate/
            ref_1/
                batch-relate-0.parquet
                report-relate.json
            ref_2/
                batch-relate-0.parquet
                report-relate.json
        mask/
            ref_1/
                batch-mask-0.csv.gz
                report-mask.json
            ref_2/
                batch-mask-0.csv.gz
                report-mask.json
    sample_2/
        align/
            ref_1.bam
            ref_2.bam
        relate/
            ref_1/
                batch-relate-0.parquet
                report-relate.json
            ref_2/
                batch-relate-0.parquet
                report-relate.json
        mask/
            ref_1/
                batch-mask-0.csv.gz
                report-mask.json
            ref_2/
                batch-mask-0.csv.gz
                report-mask.json

List of files
------------------------------------------------------------------------

The simplest manner is to list every input file explicitly. For example,
two report files could be processed with ``table`` command like this::

    seismic table sample_1/relate/ref_1/report-relate.json sample_2/mask/ref_2/report-relate.json

Glob patterns
------------------------------------------------------------------------

Listing many input files explicitly would be tedious. `Glob patterns`_
use wildcard characters to match many paths with a single expression.
This method is especially useful for matching files that have the same
names and are located in different directories. For example, to process
all report files for the reference ``ref_2``::



.. _glob patterns: https://en.wikipedia.org/wiki/Glob_(programming)
