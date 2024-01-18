
Report Formats
========================================================================

A report file is produced by the steps ``align``, ``relate``, ``pool``,
``mask``, ``cluster``, ``join``, and ``fold`` and serves three purposes:

- Record the parameters with which the step was run.
- Summarize the results of the step.
- Assist with loading the data output by the step (except ``align``).

Every report file is saved in `JSON format`_, as this format easy to
read and widely supported by software.


.. toctree::
    :maxdepth: 2

    align
    relate
    pool
    mask
    cluster
    join
    fold

.. _JSON format: https://en.wikipedia.org/wiki/JSON
