
.. _plots_get_started:

Get started
=================

Import your data into a study object
------------------------------------

.. code::

    from dreem.draw import study
    import json

    # Define the list of samples you want to load
    my_samples = ['my_sample_1', 'my_sample_2', 'my_sample_3']

    # Load the data as a list of dictionaries
    data = []
    for sample in my_samples:
        with open('{}.json'.format(sample), 'r') as f:
            data.append(json.load(f))

    # Create a study object
    my_study = study.Study(data)

    # Print the list of available samples
    print(my_study.get_samples())


Make a plot
-----------

.. code::

    study.mutation_fraction(
        sample='my_sample_1',
        reference = 'my_reference_1',
        section = 'my_section_1',
        cluster = 'pop_avg'
    )

.. note::

    We regularly update the list of available plots. Make sure that you have the last version of DREEM.
