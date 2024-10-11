
Profile: Bar graph of relationships(s) per position
--------------------------------------------------------------------------------

Profile: Input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can give any number of positional table files:

- ``relate-per-pos.csv``
- ``mask-per-pos.csv``
- ``clust-per-pos.csv``

See :ref:`graph_inputs` for more information on input files.

Profile: Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`graph_data` for more information on graph options.

Profile: Output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`graph_outputs` for more information on output files.

Profile: Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assume your output files are in the directory ``out``.

Graph the fraction of reads mutated (the mutation rate) at each position
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Command::

    seismic graph profile out

Outputs:

    - ``out/{sample}/graph/{ref}/{sect}/profile_masked_m-ratio-q0.csv`` ::

        Position,Base,Mutated
        196,A,0.03212678813091909
        197,C,0.003393658868599399
        198,A,0.010770982726367673
        199,G,
        200,T,
        ...

    - ``out/{sample}/graph/{ref}/{sect}/profile_masked_m-ratio-q0.html``

        .. image::
            profile_m-ratio.png

Graph the fraction each type of mutation at each position
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Command::

    seismic graph profile -r acgtdi out

Outputs:

    - ``out/{sample}/graph/{ref}/{sect}/profile_masked_acgtdi-ratio-q0.csv`` ::

        Position,Base,Subbed-A,Subbed-C,Subbed-G,Subbed-T,Deleted,Inserted
        196,A,0.0,0.018827932574530207,0.0016602828354697209,0.011638240863466542,0.0,0.0
        197,C,0.0012500543750758389,0.0,0.00042715097170304697,0.0017164535218205134,0.0,0.0
        198,A,0.0,0.002887023310345601,0.0005903400048243549,0.007293619411197718,0.0,0.0
        199,G,,,,,,
        200,T,,,,,,
        ...

    - ``out/{sample}/graph/{ref}/{sect}/profile_masked_acgtdi-ratio-q0.html``

        .. image::
            profile_acgtdi-ratio.png

Graph the number of informative reads at each position
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Command::

    seismic graph profile -r n --use-count out

Outputs:

    - ``out/{sample}/graph/{ref}/{sect}/profile_masked_n-count.csv`` ::

        Position,Base,Informed
        196,A,301334.2
        197,C,305746.7
        198,A,303045.7
        199,G,
        200,T,
        ...

    - ``out/{sample}/graph/{ref}/{sect}/profile_masked_n-count.html``

        .. image::
            profile_n-count.png
