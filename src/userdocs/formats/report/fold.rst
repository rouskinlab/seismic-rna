
Fold Report
--------------------------------------------------------------------------------

Fold Report: Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======================================================================================= =========
Name                                                                                    Data Type
======================================================================================= =========
Name of Sample                                                                          str
Name of Reference                                                                       str
Name of Section                                                                         str
Name of Profile                                                                         str
Quantile for normalizing ratios; must be in [0, 1]                    float
Temperature at which to predict structures (in Kelvin)                                  float
Maximum distance between two paired bases in predicted structures (use 0 for no limit)  int
Predict only the minimum free energy (MFE) structure, without any suboptimal structures bool
Maximum number of predicted structures                                                  int
Maximum % difference in energy between suboptimal structures                            float
Branches                                                                                list[str]
Time Began                                                                              str
Time Ended                                                                              str
Time Taken (minutes)                                                                    float
Version of SEISMIC-RNA                                                                  str
======================================================================================= =========

Fold Report: Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    {
        "Name of Sample": "ldi",
        "Name of Reference": "sars2_1799",
        "Name of Section": "third",
        "Name of Profile": "fse__average",
        "Quantile for normalizing ratios; must be in [0, 1]": 0.95,
        "Temperature at which to predict structures (in Kelvin)": 310.15,
        "Maximum distance between two paired bases in predicted structures (use 0 for no limit)": 0,
        "Predict only the minimum free energy (MFE) structure, without any suboptimal structures": false,
        "Maximum number of predicted structures": 20,
        "Maximum % difference in energy between suboptimal structures": 20.0,
        "Branches": [],
        "Time Began": "2023-12-29 at 22:10:46",
        "Time Ended": "2023-12-29 at 22:10:53",
        "Time Taken (minutes)": 0.12,
        "Version of SEISMIC-RNA": "0.10.0"
    }
