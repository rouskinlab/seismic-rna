

Clustering report 
+++++++++++++++++

This file is outputted by the clustering algorithm. 
It associates each read with a likelihood to belong to a cluster.

The file is a json with the following structure:


.. code:: python


    {
        # 2 clusters
        'K2_1':
        {
            'read_1': 0.5,
            'read_2': 0.8,
        },
        'K2_2':
        {
            'read_1': 0.5,
            'read_2': 0.2,
        },

        # 3 clusters
        'K3_1':
        {
            'read_1': 0.2,
            'read_2': 0.5,
        },

        'K3_2':
        {
            'read_1': 0.2,
            'read_2': 0.5,
        },

        'K3_3':
        {
            'read_1': 0.6,
            'read_2': 0.0,
        },
    }


.. note::

    Why not a csv like this?

    ========= ====== ====== ====== ====== ====== 
     read_id   K2_1   K2_2   K3_1   K3_2   K3_3  
    ========= ====== ====== ====== ====== ====== 
     read_1    0.5    0.5    0.2    0.2    0.6   
     read_2    0.8    0.2    0.5    0.5    0.0   
    ========= ====== ====== ====== ====== ====== 

