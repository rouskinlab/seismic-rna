********************************************************************************
Install
********************************************************************************


System Requirements
================================================================================

SEISMIC-RNA runs on Linux and macOS.
If you use Windows, we recommend installing and running SEISMIC-RNA using the
`Windows Subsystem for Linux (WSL)`_.


TL;DR
================================================================================

If you already have Conda_ or Mamba_ installed and know how to use it, then
type these three commands into a terminal (substitute ``mamba`` for
``conda`` if you use Mamba)::

    conda create -n seismic python=3.13
    conda activate seismic
    conda install -c conda-forge -c bioconda seismic-rna

After you have installed SEISMIC-RNA, :ref:`set_datapath`.

If this fails, or if you don't have or know how to use Conda or Mamba, then
read on.


.. _install_with_conda:

Install with Conda or Mamba
================================================================================

This section is for anyone who does not already have Conda or Mamba
installed, or is not yet comfortable using either.
It explains what Conda and Mamba are, what a virtual environment is, and
then walks you step by step through installing Conda (or Mamba), creating a
virtual environment, and installing SEISMIC-RNA into it.

Step 1: Install Conda or Mamba
--------------------------------------------------------------------------------

Conda_ is a free, popular package manager: a program that automates
downloading and installing other software, including whichever versions of
its dependencies are compatible with each other.
Mamba_ is a newer replacement for Conda, reimplemented in C++, that can
install packages much faster; the two behave the same way, so wherever this
page shows a ``conda`` command, you may instead type ``mamba``.

If you don't already have Conda or Mamba, then we recommend installing a
small version: Miniconda_ for Conda, or Miniforge for Mamba (see Mamba_ for
installation instructions).
When the installer asks if you want to initialize Conda (or Mamba), choose
yes.
Otherwise, you can initialize it later by typing the path to your ``conda``
(or ``mamba``) executable followed by ``init``.
If you install in the default location, then the command to initialize Conda
is ::

    ~/miniconda3/bin/conda init

and to initialize Mamba is ::

    ~/miniforge3/bin/mamba init

Step 2: Create a virtual environment for SEISMIC-RNA
--------------------------------------------------------------------------------

A virtual environment is an isolated space on your computer into which you
can install software -- such as SEISMIC-RNA and the specific versions of its
dependencies -- without it conflicting with any other software already on
your computer.
Create one for SEISMIC-RNA with Conda::

    conda create -n seismic python=3.13

or with Mamba::

    mamba create -n seismic python=3.13

You must indicate which version of Python to use; SEISMIC-RNA supports only
Python 3.13, so specify that version.
You can name your environment whatever you like using the ``-n`` option; in
this example, it is named ``seismic``.

.. note::

    We recommend giving your environment a short name because you will need
    to type its name every time before using it.

Step 3: Activate the virtual environment for SEISMIC-RNA
--------------------------------------------------------------------------------

Before you install SEISMIC-RNA into the virtual environment, you must
"activate" the environment using the name you gave it (which was ``seismic`` in
this tutorial).
With Conda::

    conda activate seismic

With Mamba::

    mamba activate seismic

.. warning::

    Make sure to activate the environment for SEISMIC-RNA before installing
    any packages for SEISMIC-RNA.
    If you don't, then you will instead install the packages into whichever
    environment was already active, which would not only unintentionally
    alter this environment but also fail to install the packages into the
    ``seismic`` environment.

Step 4: Install SEISMIC-RNA and its dependencies
--------------------------------------------------------------------------------

Run this command to install SEISMIC-RNA and all other software it requires using
Conda::

    conda install -c conda-forge -c bioconda seismic-rna

Or using Mamba::

    mamba install -c conda-forge -c bioconda seismic-rna

.. note::

    Conda or Mamba may fail to install SEISMIC-RNA if some of its
    dependencies are not compatible with your hardware and/or operating
    system, or if you need a version of SEISMIC-RNA that is newer than the
    one available through Conda (and Mamba).
    If that happens, then see :ref:`install_without_conda`, which explains
    how to still use Conda or Mamba for the environment and dependencies
    while installing SEISMIC-RNA itself a different way.

After you have installed SEISMIC-RNA, :ref:`set_datapath`.


.. _install_without_conda:

Install without Conda or Mamba
================================================================================

"Without Conda or Mamba" here means installing SEISMIC-RNA itself without Conda
or Mamba -- for example, if the version you need is not yet available through
either (see the warning at the top of this page).
You can still use Conda or Mamba to create your virtual environment and
install the non-Python dependencies; this section explains how to do so, as
well as how to do so without either, before explaining how to install
SEISMIC-RNA itself.

Step 1: Create and activate a virtual environment
--------------------------------------------------------------------------------

A virtual environment keeps SEISMIC-RNA and its dependencies isolated from
any other software on your computer.
You can create one with Conda_ or Mamba_::

    conda create -n seismic python=3.13
    conda activate seismic

or, if you would rather not use Conda or Mamba, with Python's built-in venv_
module (this requires that Python 3.13 already be installed; see Python_)::

    python3.13 -m venv seismic
    source seismic/bin/activate

Step 2: Install the non-Python dependencies
--------------------------------------------------------------------------------

SEISMIC-RNA depends on several pieces of non-Python software that cannot be
installed with Python's package manager ``pip``.
The easiest way to install them is with Conda or Mamba, using the following
commands:

.. note::
    Strictly speaking, you only need to install the dependencies for the
    commands you plan to use (see the "Used in" column below).
    However, we recommend installing all of them if you can, so that you do
    not run into errors caused by a missing dependency later on, after you
    have already started using SEISMIC-RNA.

=============  ==============================================================  ============================================================================================
Dependency     Command to install (``conda`` can be replaced with ``mamba``)   Used in
=============  ==============================================================  ============================================================================================
Bowtie2_       ``conda install -c conda-forge -c bioconda bowtie2``            ``align``
Fastp_         ``conda install -c conda-forge -c bioconda fastp``              ``align``
Numba_         ``conda install -c conda-forge numba>=0.67``                    all steps (see note)
RNAstructure_  ``conda install -c conda-forge -c bioconda rnastructure>=6.6``  ``fold``, ``sim fold``, ``sim total`` (default)
Samtools_      ``conda install -c conda-forge -c bioconda samtools``           ``align``, ``idmut``
seqkit_        ``conda install -c conda-forge -c bioconda seqkit>=2.13.0``     ``demult``
ViennaRNA_     ``conda install -c conda-forge -c bioconda viennarna>=2.7.2``   ``fold``, ``sim fold``, ``sim total`` (optional, but required by ``fold`` for duplex tables)
=============  ==============================================================  ============================================================================================

.. note::
    Numba_, unlike the other dependencies in this table, is a Python package
    that pip normally installs automatically along with SEISMIC-RNA; you only
    need to install it separately here if you are on a Mac with an Intel
    processor (see the :ref:`note about Intel Macs <intel_mac_numba>` below).

If you would rather not use Conda or Mamba for a dependency, or it fails to
install, then follow its link above to find instructions for installing it
manually.

.. note::
    If you install software on macOS without using Conda or Mamba, then you
    will need to manually approve each piece of software before you can run
    it.
    (This limitation is a safety feature of macOS intended to reduce the risk
    of running malware accidentally.)
    To approve the software, type ``which [program]`` (replacing
    ``[program]`` with an item from the list below) to find the path of the
    executable.
    Then in Finder, open the directory that contains the executable and
    approve it there.

Confirm that each dependency is installed by running each of these commands,
one at a time::

    which bowtie2
    which fastp
    which ct2dot  # ct2dot is part of RNAstructure
    which samtools
    which seqkit
    which RNAfold  # RNAfold is part of ViennaRNA

If the dependency is installed, then it should print out the path to it.
If it says something like ``not found``, then the dependency is not
installed.

Step 3: Install SEISMIC-RNA itself
--------------------------------------------------------------------------------

With your dependencies installed and your virtual environment activated,
install SEISMIC-RNA itself with pip, either from the Python Package Index or
from GitHub.

Option A: Install from the Python Package Index
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If you are not using Conda or Mamba to install SEISMIC-RNA, the next-best
option we recommend is to install SEISMIC-RNA from the Python Package Index
(PyPI_), which will install the latest stable release version.
In a terminal, type this command to install SEISMIC-RNA and all its Python
dependencies::

    pip install seismic-rna

.. _intel_mac_numba:

.. note::
    **Macs with Intel processors or M-series processors running an Intel
    (x86_64) build of Python via Rosetta.**
    On these computers, ``pip install seismic-rna`` fails while building a
    dependency called llvmlite, with an error that mentions ``CMake`` and
    ``Could not find a package configuration file provided by "LLVM"``.

    The reason is that Numba and llvmlite (which SEISMIC-RNA uses to speed up
    its calculations) no longer distribute versions built for Intel Macs on
    the Python Package Index, so pip tries to build llvmlite from its source
    code, which fails unless you have installed LLVM yourself.
    Macs with Apple Silicon processors running a native ARM64 build of
    Python, as well as Linux computers, are not affected.

    To fix this, install Numba with Conda or Mamba before installing
    SEISMIC-RNA::

        conda install -c conda-forge "numba>=0.67"
        pip install seismic-rna

Option B: Install from Source Code
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. warning::
    We do *not* recommend installing from source unless you need the latest
    source code (most users do not), which may be unstable or contain
    significant bugs.

In a terminal, navigate to the directory into which to install SEISMIC-RNA.
If Git_ is installed on your computer, you can download the latest source
code from the SEISMIC-RNA repository on GitHub::

    git clone https://github.com/rouskinlab/seismic-rna.git

Otherwise, open ``https://github.com/rouskinlab/seismic-rna`` in a web
browser, click "Code" then "Download ZIP", unzip the file after it has
downloaded, and move it to the directory where you want to keep the source
code.

To install SEISMIC-RNA, type ``pip install`` followed by the path of the
source code directory that you downloaded, e.g. ::

    pip install ~/Downloads/seismic-rna

If you want to be able to modify the source code after you install
SEISMIC-RNA and have those changes come into effect, then add the flag
``-e``, e.g. ::

    pip install -e ~/Downloads/seismic-rna

Otherwise, you may delete the source code after installation to save space.

After you have installed SEISMIC-RNA, :ref:`set_datapath`.


.. _set_datapath:

Set the DATAPATH Environment Variable
================================================================================

RNAstructure_ requires an environment variable called ``DATAPATH`` to point to
the directory of thermodynamic data tables.
See https://rna.urmc.rochester.edu/Text/Thermodynamics.html for details.
SEISMIC-RNA should be able to guess the correct ``DATAPATH`` if RNAstructure was
installed manually from the website or with Conda or Mamba, but it will log a
warning message to inform you that it had to guess.
To suppress this warning, you can create an environment variable called
``DATAPATH`` on your system.
To find the location of the data tables for RNAstructure, type ::

    seismic -q datapath

This command should print a message that begins with ``DATAPATH=``.
Add this entire line (including ``DATAPATH=``) to the end of your shell RC file:
``~/.bashrc`` on most Linux systems, ``~/.zshrc`` on most macOS systems.
Restart your terminal for the changes to take effect.
After restarting the terminal, confirm ``DATAPATH`` is set by typing ::

    echo $DATAPATH

which should print out the path to the data tables that you found previously.
Now the ``DATAPATH`` will be set automatically every time you open the terminal,
unless you remove or edit that line in your shell RC file.


.. _install_update:

Update to Another Version of SEISMIC-RNA
================================================================================

If you have already installed SEISMIC-RNA, follow these steps to install a
different version.

.. note::

    If the version you are updating to has substantially different
    dependencies from the one you have installed -- for example, if it
    requires a newer version of Python -- then updating may fail.
    If that happens, create a new virtual environment and follow the
    instructions above (:ref:`install_with_conda` or
    :ref:`install_without_conda`) as if installing SEISMIC-RNA for the first
    time.

Update SEISMIC-RNA to the latest stable version
--------------------------------------------------------------------------------

Type this to install with Conda (if the version you want to install has been
released on Bioconda)::

    conda update -c conda-forge -c bioconda seismic-rna

or with Mamba::

    mamba update -c conda-forge -c bioconda seismic-rna

Type this to install with ``pip`` (if you don't want to use Conda/Mamba or if
the version you want to install has not been released on Bioconda)::

    pip install -U seismic-rna

Install a specific version of SEISMIC-RNA
--------------------------------------------------------------------------------

Every version of SEISMIC-RNA has three parts -- x.y.z -- where x is the major
version, y is the minor version, and z is the patch (this system is known as
semantic versioning).

Type this to install with Conda (if the version you want to install has been
released on Bioconda)::

    conda install -c conda-forge -c bioconda seismic-rna=x.y.z

or with Mamba::

    mamba install -c conda-forge -c bioconda seismic-rna=x.y.z

Type this to install with ``pip`` (if you don't want to use Conda/Mamba or if
the version you want to install has not been released on Bioconda)::

    pip install -U seismic-rna==x.y.z

.. note::

    When specifying the version, use ``=`` with Conda or Mamba and ``==``
    with pip.


.. _test_seismicrna:

Test SEISMIC-RNA
================================================================================

SEISMIC-RNA comes with hundreds of tests to verify that it is working properly
on your system.

.. note::
    Running the full test suite can take more than 30 minutes, so it is
    entirely optional. Every release of SEISMIC-RNA is already tested to
    confirm that all tests pass before it is published, so you do not need
    to run the tests yourself unless you suspect a problem with your own
    installation.

Step 1: Run SEISMIC-RNA's testing suite
--------------------------------------------------------------------------------

To run all the tests, type this::

    seismic test

To monitor the tests as they run, you can use verbose mode (option ``-v``).
In verbose mode, as each test finishes, it will print ``.`` if it succeeds,
``F`` if it fails, ``E`` if it errs, and ``s`` if it was skipped::

    seismic test -v

To print out the name of each test as it runs and check which tests succeed and
fail, you can use double-verbose mode::

    seismic test -vv

Step 2: Interpret the test results
--------------------------------------------------------------------------------

Regardless of the verbosity, if all tests succeed, then it will print a message
similar to this::

    Ran 1651 tests in 1176.813s

    OK

Otherwise, it will print the number of tests that failed and a message about
each failure.
If this happens, then first follow :ref:`install_update` to ensure you are using
the latest version of SEISMIC-RNA and its dependencies.
If your problem persists, then please report an issue (see :doc:`../issues` for
instructions).


.. _Conda: https://docs.conda.io/en/latest/
.. _Mamba: https://mamba.readthedocs.io/en/latest/
.. _Bioconda: https://bioconda.github.io/
.. _Git: https://git-scm.com/
.. _Miniconda: https://docs.anaconda.com/miniconda/
.. _pip: https://pip.pypa.io/en/stable/
.. _Python: https://www.python.org/downloads/
.. _venv: https://docs.python.org/3/library/venv.html
.. _Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _Fastp: https://github.com/OpenGene/fastp
.. _Numba: https://numba.pydata.org/
.. _RNAstructure: https://rna.urmc.rochester.edu/RNAstructure.html
.. _Samtools: https://www.htslib.org/
.. _seqkit: https://bioinf.shenwei.me/seqkit/
.. _ViennaRNA: https://www.tbi.univie.ac.at/RNA/
.. _PyPI: https://pypi.org/project/seismic-rna/
.. _Anaconda: https://anaconda.org/bioconda/seismic-rna
.. _Windows Subsystem for Linux (WSL): https://learn.microsoft.com/en-us/windows/wsl
