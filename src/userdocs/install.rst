********************************************************************************
Install
********************************************************************************


Option 1: Quick installation (if you already have Python_, pip_, and Conda_)
================================================================================

Run these commands in a terminal::

    # Create and activate a new Conda environment.
    conda create -n seismic python=3.12
    conda activate seismic
    # Install the non-Python dependencies.
    conda install -y -c bioconda -c conda-forge bowtie2 fastqc rnastructure samtools
    # Install SEISMIC-RNA and its Python dependencies.
    pip install seismic-rna

You can name your environment whatever you like using the ``-n`` option; in this
example, it is named ``seismic``.
SEISMIC-RNA is compatible with Python version 3.10 and later.


Option 2: Installation with Conda
================================================================================

Install Conda and pip
--------------------------------------------------------------------------------

We highly recommend installing SEISMIC-RNA into a virtual environment to spare
yourself future frustration.
Conda_ is a popular tool for managing virtual environments, especially (but not
exclusively) for Python-based software.
We recommend using the Miniconda_ installer, which installs both Conda and pip.
When the installer asks if you want to initialize Conda, choose yes.
If you do not, you can initialize Conda later by typing the path to your Conda
executable followed by ``init``, e.g. ::

    ~/miniconda3/bin/conda init

Create a Conda environment for SEISMIC-RNA
--------------------------------------------------------------------------------

Once Conda and pip are installed, create a new virtual environment into which
SEISMIC-RNA and all other necessary software will go::

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

Install the non-Python dependencies
--------------------------------------------------------------------------------

SEISMIC-RNA requires several other pieces of software:

============ =================================================================================================
Software     SEISMIC-RNA commands that use the software
============ =================================================================================================
Bowtie2      ``seismic align``; ``seismic wf``
FastQC       ``seismic align`` (without ``--no-fastqc``); ``seismic wf`` (without ``--no-fastqc``)
RNAstructure ``seismic fold``; ``seismic wf`` (with ``--fold``); ``seismic +sim fold``; ``seismic +sim total``
Samtools     ``seismic align``; ``seismic relate``; ``seismic wf``
============ =================================================================================================

Install these pieces of software with this command::

    conda install -y -c bioconda -c conda-forge bowtie2 fastqc rnastructure samtools

Install SEISMIC-RNA and its Python dependencies
--------------------------------------------------------------------------------

Finally, install SEISMIC-RNA and all of the Python packages that it requires::

    pip install seismic-rna


Option 3: Installation without Conda
================================================================================

Although we highly recommend using Conda_ or other software that supports virual
environments, you can also install SEISMIC-RNA without it.

Install Python and pip
--------------------------------------------------------------------------------

First, if Python_ is not installed, then install the latest version.
Confirm that Python version 3.10 or later and pip_ are installed::

    python3 --version
    pip3 --version

Install SEISMIC-RNA and its Python dependencies
--------------------------------------------------------------------------------

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
and have those changes come into effect, then add the flag ``-e``.
Otherwise, to save space, you may delete the source code after installation.


Option 4: Upgrading (if you already have SEISMIC-RNA installed)
================================================================================

To upgrade SEISMIC-RNA to the latest version, type ::

    pip install -U seismic-rna


To install a specific version ``x.y.z``, type ::

    pip install seismic-rna==x.y.z


Troubleshooting installation
================================================================================


.. _conda: https://docs.conda.io/en/latest/
.. _git: https://git-scm.com/
.. _miniconda: https://docs.anaconda.com/miniconda/
.. _pip: https://pip.pypa.io/en/stable/
.. _python: https://www.python.org/downloads/
