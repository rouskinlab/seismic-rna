********************************************************************************
Install
********************************************************************************


Option 1: Quick installation (if you already have Conda_)
================================================================================

SEISMIC-RNA is available from the Bioconda_ channel and compatible with Python
3.10 and later.
To install, type this into a terminal::

    conda create -n seismic python=3.12
    conda activate seismic
    conda install -c bioconda -c conda-forge seismic-rna

If installing with Conda fails, then try :ref:`install_without_conda`.


Option 2: Step-by-step installation with Conda
================================================================================

Install Conda
--------------------------------------------------------------------------------

Installing SEISMIC-RNA is easiest with Conda_, a popular package manager.
If you don't already have Conda, then we recommend installing a small version of
it called Miniconda_.
When the installer asks if you want to initialize Conda, choose yes.
Otherwise, you can initialize Conda later by typing the path to your ``conda``
executable followed by ``init``, e.g. ::

    ~/miniconda3/bin/conda init

Create a Conda environment for SEISMIC-RNA
--------------------------------------------------------------------------------

Once Conda is installed, create a new virtual environment into which SEISMIC-RNA
and all other necessary software will go::

    conda create -n seismic python=3.12

You can name your environment whatever you like using the ``-n`` option; in this
example, it is named ``seismic``.

.. note::

    We recommend giving your environment a short name because you will need to
    type its name every time before using it.

You must indicate which version of Python to use; we recommend the most recent
stable release (currently version 3.12), though SEISMIC-RNA is compatible with
version 3.10 and later.

Activate the Conda environment for SEISMIC-RNA
--------------------------------------------------------------------------------

Before you install SEISMIC-RNA into the Conda environment, you must "activate"
the environment by typing ``conda activate`` followed by its name, e.g. ::

    conda activate seismic

.. warning::

    Make sure to activate the environment for SEISMIC-RNA before installing any
    packages for SEISMIC-RNA.
    If you don't, then you will instead install the packages into whichever
    environment was already active, which would not only unintentionally alter
    this environment but also fail to install the packages into the ``seismic``
    environment.










ADD A SEGMENT ON WHAT TO DO IF CONDA INSTALL FAILS BUT YOU CAN STILL USE CONDA







Install SEISMIC-RNA and its dependencies
--------------------------------------------------------------------------------

Run this command to install SEISMIC-RNA and all other software it requires::

    conda install -c bioconda -c conda-forge seismic-rna

If installing with Conda fails, then try :ref:`install_without_conda`.


.. _install_without_conda:

Option 3: Step-by-step installation without Conda
================================================================================

Although Conda is the easiest means to install SEISMIC-RNA, it is not necessary.

Install Python and pip
--------------------------------------------------------------------------------

First, if Python_ is not installed, then install the latest version.
Confirm that Python version 3.10 or later and pip_ are installed::

    python --version
    pip --version

Install SEISMIC-RNA and its Python dependencies
--------------------------------------------------------------------------------

Option 3a: Install SEISMIC-RNA from the Python Package Index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend installing SEISMIC-RNA from the Python Package Index, which will
download the latest stable version that has been released.
In a terminal, type this command to install it and all its Python dependencies::

    pip install seismic-rna

Option 3b: Install SEISMIC-RNA from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We do *not* recommend installing from GitHub unless you need the latest source
code (most users do not), which may be unstable or contain significant bugs.
In a terminal, navigate to the directory into which to install SEISMIC-RNA.
If Git_ is installed on your computer, then clone the GitHub repository::

    git clone https://github.com/rouskinlab/seismic-rna.git

Otherwise, open ``https://github.com/rouskinlab/seismic-rna`` in a web browser,
click "Code" then "Download ZIP", unzip the file after it has downloaded, and
move it to the directory where you want to keep the source code.

To install SEISMIC-RNA, type ``pip install`` followed by the path of the source
code directory that you downloaded, e.g. ::

    pip install ~/Downloads/seismic-rna

If you want to be able to modify the source code after you install SEISMIC-RNA
and have those changes come into effect, then add the flag ``-e``, e.g. ::

    pip install -e ~/Downloads/seismic-rna

Otherwise, you may delete the source code after installation to save space.

Install the non-Python dependencies of SEISMIC-RNA
--------------------------------------------------------------------------------

SEISMIC-RNA requires several other pieces of software that cannot be installed
using pip:

============ =================================================================================================
Software     SEISMIC-RNA commands that use the software
============ =================================================================================================
Bowtie2      ``seismic align``; ``seismic wf``
FastQC       ``seismic align`` (without ``--no-fastqc``); ``seismic wf`` (without ``--no-fastqc``)
RNAstructure ``seismic fold``; ``seismic wf`` (with ``--fold``); ``seismic +sim fold``; ``seismic +sim total``
Samtools     ``seismic align``; ``seismic relate``; ``seismic wf``
============ =================================================================================================

You only need to install the software for the steps you will use: for example,
if you won't model RNA structures (using ``seismic fold``, ``seismic wf`` with
the ``--fold`` option, ``seismic +sim fold``, or ``seismic +sim total``), then
you don't need to install RNAstructure.
Although all of these packages are available through Conda, if you are not using
Conda (or if the Conda installation fails) then you can also install them using
the installation instructions in their documentation:

- Bowtie2_ (make sure to install Bowtie2, not Bowtie)
- FastQC_
- RNAstructure_
- Samtools_


.. _install_update:

Option 4: Update to another version (if you already installed SEISMIC-RNA)
================================================================================

Update SEISMIC-RNA to the latest stable version
--------------------------------------------------------------------------------

Type this if you had initially installed SEISMIC-RNA with Conda::

    conda update -c bioconda -c conda-forge seismic-rna

or this if you had initially installed it with pip::

    pip install -U seismic-rna


Install a specific version of SEISMIC-RNA
--------------------------------------------------------------------------------

Type this if you had initially installed SEISMIC-RNA with Conda::

    conda install -c bioconda -c conda-forge seismic-rna=x.y.z

or this if you had initially installed it with pip::

    pip install seismic-rna==x.y.z

.. note::

    When specifying the version, use ``=`` with Conda and ``==`` with pip.

.. _install_update_depend:

Update the dependencies of SEISMIC-RNA
--------------------------------------------------------------------------------

Type this if you had initially installed SEISMIC-RNA with Conda::

    conda update -c bioconda -c conda-forge --all

If you had initially installed it with pip, then first check that you have the
latest versions of Bowtie2_, FastQC_, RNAstructure_, and Samtools_ by typing
each of these commands::

    bowtie2 --version
    fastqc --version
    Fold --version
    samtools --version

and install the most recent versions as necessary.
Then update the dependencies of SEISMIC-RNA to the latest compatible versions::

    pip install -U seismic-rna

.. warning::

    If you update a package with pip, then it will install the latest version
    *without* checking if the package is a dependency of any other package.
    For example, SEISMIC-RNA depends on NumPy but is incompatible with the
    latest version of NumPy; ``pip install -U numpy`` will install the latest
    version of NumPy, regardless, which will cause SEISMIC-RNA to crash.

    However, if you update a package, then pip will also update its dependencies
    and *will* ensure that the versions of the dependencies it installs are
    compatible with the main package being updated.
    Thus, the safe way to update all dependencies of a package using pip is to
    update a package that is *not* a dependency of another package.
    Assuming you have not installed any packages that depend on SEISMIC-RNA,
    you can safely type ``pip install -U seismic-rna``, which will install the
    latest version of SEISMIC-RNA and the latest *compatible* version of each of
    its dependencies.


Test SEISMIC-RNA after installing
================================================================================

SEISMIC-RNA comes with hundreds of tests to verify that it is working properly
on your system.
We recommend running the tests each time you install or update SEISMIC-RNA.

Run SEISMIC-RNA's testing suite
--------------------------------------------------------------------------------

To run all the tests, type this::

    seismic +test

To monitor the tests as they run, you can use verbose mode (option ``-v``).
In verbose mode, as each test finishes, it will print ``.`` if it succeeds,
``F`` if it fails, ``E`` if it errs, and ``s`` if it was skipped::

    seismic +test -v

To print out the name of each test as it runs and check which tests succeed and
fail, you can use double-verbose mode::

    seismic +test -vv

Interpret the test results
--------------------------------------------------------------------------------

Regardless of the verbosity, if all tests succeed, then it will print a message
similar to this::

    Ran 903 tests in 196.699s

    OK

Otherwise, it will print the number of tests that failed and a message about
each failure.
If this happens, then first check that you are using the latest version of
SEISMIC-RNA by typing ::

    seismic --version

and checking what the version is on PyPI_ or Anaconda_.
If your installed version is not the latest, then try :ref:`install_update`.
If you already have the latest version, then the problem could be caused by one
of the dependencies, so try to :ref:`install_update_depend`.
If that does not solve your problem either, then please report a bug (see
:doc:`./issues` for instructions).


.. _Conda: https://docs.conda.io/en/latest/
.. _Bioconda: https://bioconda.github.io/
.. _Git: https://git-scm.com/
.. _Miniconda: https://docs.anaconda.com/miniconda/
.. _pip: https://pip.pypa.io/en/stable/
.. _Python: https://www.python.org/downloads/
.. _Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _RNAstructure: https://rna.urmc.rochester.edu/RNAstructure.html
.. _Samtools: https://www.htslib.org/
.. _PyPI: https://pypi.org/project/seismic-rna/
.. _Anaconda: https://anaconda.org/bioconda/seismic-rna
