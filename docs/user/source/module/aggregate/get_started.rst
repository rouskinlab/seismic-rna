
Get started
++++++++++++++++


Let's aggregate ``sample_1``.

You need:
    - the path to the bitvectors of ``sample_1`` stored as shown in the :ref:`aggregate_io` section.
    - the reference genome under the `fasta format <https://en.wikipedia.org/wiki/FASTA_format>`_ ``reference.fasta``.

You get:
    -  ``sample_1.json`` under the :ref:`dreem_output` format.



CLI
---------

Run the command:

::
    
    dreem aggregate --bv-files path/to/sample_1 —-fasta path/to/reference.fasta 



Python
------------

Run the following code:

:: 

    from dreem import aggregate
    aggregate.run(
        bv_files =('path/to/sample_1',), 
        fasta='path/to/reference.fasta'
        )


See the :ref:`aggregate_reference` for more details.
