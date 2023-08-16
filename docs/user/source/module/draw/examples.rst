
Examples
++++++++

.. code:: bash

    # Plots mutation_fraction, mutation_fraction_identity, and base_coverage
    dreem draw --inpt my_test_sample.json 
    
    # Plots only the mutation_fraction plot
    dreem draw --inpt my_test_sample.json --mutation_fraction 

    # Plots only the specified regions for each reference
    dreem draw --inpt my_test_sample.json -c reference_1 17 69 -c reference_2 18 21

    # Get a structured output 
    dreem draw --inpt my_test_sample.json --flat 0 

    # Give name to the sections in the output names
    dreem draw --inpt my_test_sample.json --flat 0 --library library.csv