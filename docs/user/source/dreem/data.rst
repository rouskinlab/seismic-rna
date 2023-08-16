Download a test dataset
-----------------------

DREEM pipeline dataset
***********************

The DREEM input dataset is available at [insert link here].

DREEM.draw dataset
********************

Use a test dataset to play with the visualization tools.

.. code:: python

    import dreem
    study = dreem.draw.Study()
    study.df = dreem.draw.load_dataset()
