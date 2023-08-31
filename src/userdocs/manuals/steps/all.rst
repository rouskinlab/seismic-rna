
All steps with one command
------------------------------------------------------------------------

.. note::
    ``seismic all`` accepts FASTQ, BAM, relate/mask/cluster report, and
    table files and directories as inputs.

From BAM, report, and/or table file(s)::

    seismic all refs.fa out/sample/align/Ref.bam out/sample/*/report-*.json out/sample/table/*/*.csv


.. note::
    Only the align, relate, mask, and table steps run by default. Enable
    clustering by specifying ``--max-clusters`` (``-k``) followed by the
    maximum number of clusters to attempt. Enable structure prediction
    with the flag ``--fold``.
