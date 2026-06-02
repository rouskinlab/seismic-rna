********************************************************************************
seismic test
********************************************************************************


Purpose
================================================================================

``seismic test`` runs SEISMIC-RNA's built-in test suite using Python's
``unittest`` framework.
Use it after installation or after upgrading to verify that the software is
working correctly on your system.


Inputs
================================================================================

No input files are required.
The test suite is self-contained and uses synthetic data generated internally.


Outputs
================================================================================

Test results are printed to standard output.
Exit code 0 means all tests passed; a non-zero exit code means one or more
tests failed.


Quick example
================================================================================

Run the full test suite::

    seismic test

Run with a progress character (``.``/``F``/``E``/``s``) per test::

    seismic test -v

Run with one line per test (its name and result)::

    seismic test -vv


Options
================================================================================

``-v`` / ``-vv``
    Increase test output verbosity: ``-v`` prints a progress character
    (``.``/``F``/``E``/``s``) per test; ``-vv`` prints one line per test.
    (Controlled by Python ``unittest``, not SEISMIC-RNA's logging system —
    see :doc:`/use/logging`.)

The auto-generated :doc:`/cli` lists every option with its current default.


Common unexpected results
================================================================================

Tests skipped
    Some tests require optional dependencies (e.g. RNAstructure, ViennaRNA,
    RNArtistCore) and are automatically skipped if those are not installed.
    Skipped tests are not failures.

Tests fail after upgrade
    Run ``seismic migrate`` if you upgraded from 0.24 to 0.25 and are using
    existing output files.
    Otherwise, reinstall SEISMIC-RNA and its dependencies.


See also
================================================================================

- :doc:`/use/logging` — logging and verbosity flags
- :doc:`/use/utility/migrate` — migrate output directories after upgrading
