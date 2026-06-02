
Fold Report
-----------

Fold Report: Fields
^^^^^^^^^^^^^^^^^^^

=============================================================================================================================================================================================================== ==============
Name                                                                                                                                                                                                            Data Type     
=============================================================================================================================================================================================================== ==============
Sample                                                                                                                                                                                                          str           
Branches                                                                                                                                                                                                        dict[str, str]
Reference                                                                                                                                                                                                       str           
Region                                                                                                                                                                                                          str           
Profile                                                                                                                                                                                                         str           
Only generate the fold command and input files; do not run folding                                                                                                                                              bool          
Model RNA structures using RNAstructure (Fold/ShapeKnots) or ViennaRNA (RNAfold/RNAsubopt); auto selects RNAstructure for DMS and ViennaRNA for other probes                                                    str           
Use this method to incorporate reactivities into folding energies. auto selects Cordero for DMS and Eddy for other probes; Eddy requires --fold-backend=ViennaRNA; Cordero requires --fold-backend=RNAstructure str           
Slope (kcal/mol) for SHAPE reactivities using Deigan method; used only with --fold-energy-method=Deigan                                                                                                         float         
Intercept (kcal/mol) for SHAPE reactivities using Deigan method; used only with --fold-energy-method=Deigan                                                                                                     float         
Normalize and winsorize reactivities to this quantile for folding                                                                                                                                               float         
Allow isolated (non-stacked) base pairs when folding                                                                                                                                                            bool          
Predict structures at this temperature (Celsius)                                                                                                                                                                float         
Limit base pair distances to this number of bases (0 for no limit)                                                                                                                                              int           
Predict only the minimum free energy (MFE) structure                                                                                                                                                            bool          
Output at most this many structures (overriden by --fold-mfe)                                                                                                                                                   int           
Stop outputting structures when the % difference in energy exceeds this value (overriden by --fold-mfe)                                                                                                         float         
Maximum absolute energy difference (kcal/mol) from the MFE for suboptimal structures output by RNAsubopt (overriden by --fold-mfe)                                                                              float         
Predict pseudoknotted structures (requires --fold-backend=RNAstructure; uses ShapeKnots when set, Fold otherwise)                                                                                               bool          
Checksum of the ViennaRNA command file (SHA-512)                                                                                                                                                                str           
Checksum of the constraints file (SHA-512)                                                                                                                                                                      str           
Checksum of the Eddy paired-prior file (SHA-512)                                                                                                                                                                str           
Checksum of the Eddy unpaired-prior file (SHA-512)                                                                                                                                                              str           
Time began                                                                                                                                                                                                      str           
Time ended                                                                                                                                                                                                      str           
Time taken (minutes)                                                                                                                                                                                            float         
Version of SEISMIC-RNA                                                                                                                                                                                          str           
=============================================================================================================================================================================================================== ==============

Fold Report: Example
^^^^^^^^^^^^^^^^^^^^

::

    {
        "Sample": "sample1",
        "Branches": {
            "align": "",
            "idmut": "",
            "filter": "",
            "fold": ""
        },
        "Reference": "myref",
        "Region": "full",
        "Profile": "full__average",
        "Only generate the fold command and input files; do not run folding": false,
        "Model RNA structures using RNAstructure (Fold/ShapeKnots) or ViennaRNA (RNAfold/RNAsubopt); auto selects RNAstructure for DMS and ViennaRNA for other probes": "RNAstructure",
        "Use this method to incorporate reactivities into folding energies. auto selects Cordero for DMS and Eddy for other probes; Eddy requires --fold-backend=ViennaRNA; Cordero requires --fold-backend=RNAstructure": "Cordero",
        "Slope (kcal/mol) for SHAPE reactivities using Deigan method; used only with --fold-energy-method=Deigan": 1.8,
        "Intercept (kcal/mol) for SHAPE reactivities using Deigan method; used only with --fold-energy-method=Deigan": -0.6,
        "Normalize and winsorize reactivities to this quantile for folding": 0.95,
        "Allow isolated (non-stacked) base pairs when folding": false,
        "Predict structures at this temperature (Celsius)": 37.0,
        "Limit base pair distances to this number of bases (0 for no limit)": 0,
        "Predict only the minimum free energy (MFE) structure": false,
        "Output at most this many structures (overriden by --fold-mfe)": 20,
        "Stop outputting structures when the % difference in energy exceeds this value (overriden by --fold-mfe)": 20.0,
        "Maximum absolute energy difference (kcal/mol) from the MFE for suboptimal structures output by RNAsubopt (overriden by --fold-mfe)": 1.0,
        "Predict pseudoknotted structures (requires --fold-backend=RNAstructure; uses ShapeKnots when set, Fold otherwise)": false,
        "Checksum of the ViennaRNA command file (SHA-512)": "",
        "Checksum of the constraints file (SHA-512)": "",
        "Checksum of the Eddy paired-prior file (SHA-512)": "",
        "Checksum of the Eddy unpaired-prior file (SHA-512)": "",
        "Time began": "2026-05-31 at 12:59:24",
        "Time ended": "2026-05-31 at 12:59:24",
        "Time taken (minutes)": 0.01,
        "Version of SEISMIC-RNA": "0.25.3"
    }
