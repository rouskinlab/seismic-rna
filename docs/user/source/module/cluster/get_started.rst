
Get started
++++++++++++++++

Let's cluster the bitvectors of ``sample_1``.

You need:
    - the path to the bitvectors of ``sample_1`` stored as shown in the :ref:`cluster_io` section.

You get:
    -  the likelihod of each read belonging to each cluster in ``best_cluster_reads.json`` under the :ref:`clustering_output` format.



CLI
---------

Run the command:

::
    
   dreem cluster --mp-report /path/to/report.json --out-dir /path/to/output_dir


Python
------------

Run the following code:

:: 

    from dreem import cluster

    # Important for multiprocessing
    if __name__ == '__main__':
        cluster.run(
            mp_report=['/path/to/report.json'],
            out_dir='/path/to/output_dir'
            )

See the :ref:`cluster_reference` for more details.