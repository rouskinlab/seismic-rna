
Graph: Plot data from tables and/or structures and compare samples
--------------------------------------------------------------------------------

Graph: Types of graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SEISMIC-RNA can create eight types of graphs:

======== ================================================================================
Graph    Description
======== ================================================================================
profile  Bar graph of relationships(s) (e.g. mutation rate or read coverage) per position
delprof  Bar graph of differences between two profiles per position
scatter  Scatter plot comparing two profiles
corroll  Rolling correlation/comparison of two profiles
histpos  Histogram of relationship(s) per position
histread Histogram of relationship(s) per read
roc      Receiver operating characteristic (ROC) curve comparing a profile to a structure
aucroll  Rolling area under the ROC curve (AUC-ROC) comparing a profile to a structure
======== ================================================================================

See :doc:`../graph/index` for details and instructions.

Graph: Default graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can make six graphs when you run the entire workflow (``seismic wf``):

- Profile of the fraction of reads mutated (the mutation rate) at each position.
- Profile of the fraction of each type of mutation at each position.
- Profile of the number of reads with information at each position.
- Histogram of the number of positions mutated in each read.
- ROC for each predicted structure (with ``--fold``).
- Rolling AUC-ROC for each predicted structure (with ``--fold``).

These graphs are created by default in ``seismic wf``.
To make additional types of graphs, see :doc:`../graph/index`.
