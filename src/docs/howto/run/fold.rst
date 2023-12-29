
Fold: Predict RNA secondary structures using mutation rates
--------------------------------------------------------------------------------

Fold: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fold input file: Reference sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You need one file of reference sequences in FASTA format (for details on this
format, see :doc:`../../formats/data/fasta`).
If your file has characters or formatting incompatible with SEISMIC-RNA, then
you can fix it using the :doc:`../cleanfa` tool.

Fold input file: Mask or Cluster positional table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can give any number of positional table files of masked or clustered reads
(``mask-per-pos.csv`` or ``clust-per-pos.csv``, respectively) as inputs for the
Fold step.
See :doc:`../inputs` for ways to list multiple files.
(SEISMIC-RNA will not crash if you give other type of table files, such as a
``relate-per-pos.csv`` or ``mask-per-read.csv.gz`` file, but will ignore them.)

To predict structures using the mutational profiles in all valid tables in the
directory ``{out}``, you could use the command ::

    seismic fold {refs.fa} {out}

where ``{refs.fa}`` is the path to the file of reference sequences.

Fold: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fold setting: Define sections
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can predict structures of the full reference sequences or specific sections.
See :doc:`../sections` for ways to define sections.

Defining sections in ``seismic fold`` works identically to ``seismic mask`` but
accomplishes a very different purpose.
Sections in ``seismic fold`` determine for which parts of the reference sequence
to predict structures.
Sections in ``seismic mask`` determine for which parts of the reference sequence
to use mutational data.
SEISMIC-RNA allows these sections to be different.
There are several common scenarios:

- The section you are folding matches the section for which you have data.
  For example, you could have mutationally profiled a full transcript and now
  want to predict the structure of the full transcript using the data from the
  full mutational profile.
- You are folding a section that contains and is longer than the section for
  which you have data.
  For example, you could have mutationally profiled a short amplicon from a much
  longer transcript; and after clustering that amplicon, you want to model each
  alternative structure of the long transcript while using the short mutational
  profile of each cluster to guide the structure predictions.
- You are folding a short section that is contained by a longer section for
  which you have mutational profiling data.
  For example, you could have mutationally profiled a full transcript and now
  want to predict the structure of a small part of the transcript that you are
  reasonably sure does not interact with any other part of the transcript.

Fold setting: Quantile for normalization
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Folding requires that the mutation rates be normalized to the interval [0, 1].
See :doc:`../normalize` for ways to normalize mutation rates.

Fold: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All output files go into the directory ``{out}/{sample}/fold/{ref}/{sect}``,
where ``{out}`` is the output directory, ``{sample}`` is the sample, ``{ref}``
is the reference, and ``{sect}`` is the section you *folded* (**not** that from
which the data came).
The files for each predicted structure are named ``{sect}__{profile}``, where
``{sect}`` is the section from which the *data* came (**not** that which you
folded) and ``{profile}`` is the mutational profile of those data, which can be
``average`` (ensemble average) or ``cluster-{n}-{i}`` (where ``{n}`` is the
number of clusters and ``{i}`` is the cluster number).

Fold output file: Connectivity table
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The primary output is a connectivity table file.
For details on this format, see :doc:`../../formats/data/ct`.

.. _fold_db:

Fold output file: Dot-bracket structure
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Fold step also outputs the structures in dot-bracket format, which you can
copy-paste into RNA drawing software such as `VARNA`_.
For details on this format, see :doc:`../../formats/data/db`.

Fold output file: VARNA color file
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The Fold step also outputs the normalized mutation rates in VARNA color format,
which you can import into the RNA drawing software `VARNA`_.
For details on this format, see :doc:`../../formats/data/varna-color`.

Fold: Visualize structures in VARNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`VARNA`_ is a third-party application for drawing RNA structures.
To draw a structure from SEISMIC-RNA in VARNA:

1.  Install (if needed) and launch VARNA.
2.  Open your dot-bracket file (see :ref:`fold_db`) in a text editor.
3.  Right-click the drawing canvas, select "File" > "New...", and copy-paste the
    sequence and dot-bracket structure.
4.  Adjust the layout of the structure by clicking and dragging.
5.  To color the bases by their mutation rates, right-click the drawing canvas,
    select "Display" > "Color map" > "Load values...", copy-paste the path to
    your VARNA color file into the box or click "Choose file" and navigate to
    your VARNA color file, and click "OK" to load the file.
6.  To customize the colors, select "Display" > "Color map" > "Style...":

    - Drag a color bar to adjust its location.
    - Click the square below a color bar to change its color.
    - Click the X below the square to delete the color.
    - Click anywhere on the color spectrum to create a new color bar.

    We recommend setting the color for missing data (-1) to white or light gray
    and using a continuous (not discrete) color scale for the mutation data.

.. _VARNA: https://varna.lisn.upsaclay.fr/
