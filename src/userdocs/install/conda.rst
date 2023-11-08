
Installing with Conda
========================================================================

This section uses Conda to install SEISMIC-RNA in a virtual environment.
If you are not familiar with virtual environments, read this primer on
:ref:`virtual-envs`. If you do not have a working installation of Conda
on your system (or are not sure if you do), then see the section of this
tutorial about :ref:`install-prereqs` for guidance.


Create and activate your virtual environment
------------------------------------------------------------------------

Make a new virtual environment into which you will install SEISMIC-RNA.
Give it a short, memorable name, such as ``seismic``, and process Python
version 3.10, which is currently the only supported version. The rest of
this tutorial assumes that the virtual environment is named ``seismic``,
so if you have chosen a different name, then type that name instead::

    conda create -n seismic python=3.10

Activate your virtual environment with this command::

    conda activate seismic

Note that this environment **must** be active whenever you install, use,
or update SEISMIC-RNA. To deactivate it after you are finished, type::

    conda deactivate

Note that you do not need to type ``seismic`` after ``deactivate``.


Install the non-Python dependencies
------------------------------------------------------------------------

SEISMIC-RNA has four dependencies that are not Python packages. They are
all available to install using Conda, but this installation process
seems prone to failing on some systems, such as Macbook computers with
M1 and M2 processors. If you encounter errors, you can also install each
dependency by following its own installation instructions.

The four non-Python dependencies are

- `Samtools`_ ≥ v1.17 (required)
- `FastQC`_ ≥ 0.12.1 (optional: only for FASTQ quality control)
- `Bowtie2`_ ≥ v2.5.1 (optional: only for FASTQ alignment)
- `RNAstructure`_ ≥ v6.3 (optional: only for RNA structure prediction)


Installing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can install them all at once by running the following command (note
that it is long and you may need to scroll horizontally to see it all)::

    conda install -y -c bioconda -c conda-forge samtools=1.17 fastqc=0.12.1 bowtie2=2.5.1 rnastructure=6.3

Or, if you do not want to install all of them, then you can choose which
to install by running the corresponding individual command(s)::

    # Samtools
    conda install -y -c bioconda -c conda-forge samtools=1.17

    # FastQC
    conda install -y -c bioconda -c conda-forge fastqc=0.12.1

    # Bowtie2
    conda install -y -c bioconda -c conda-forge bowtie2=2.5.1

    # RNAstructure
    conda install -y -c bioconda -c conda-forge rnastructure=6.3


Verifying
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To confirm that the dependencies have been installed successfully, type
the following commands for each dependency that you want to check::

    # Samtools
    which samtools

    # FastQC
    which fastqc

    # Bowtie2
    which bowtie2

    # RNAstructure
    which ct2dot

If a dependency has been installed properly, then running its above
command will print a message that looks similar to this::

    /home/your-name/miniconda3/envs/seismic/bin/name-of-dependency

If it prints an error message that the depdendency was not found, then
the program has not been installed.


Troubleshooting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you encounter errors during installation with Conda, then first
double check that you have entered the command exactly, especially that

- the command includes ``-c bioconda`` and ``-c conda-forge``
- the version of the program in the command is correct
- there are no extra spaces or missing spaces
- every word is spelled correctly

If you had mistyped the command, then just re-run the command with the
correct spelling. If you have verified that the command is correct and
it is still not installing, then try to install the dependency using the
instructions on its website (instead of with Conda):

- Samtools: https://www.htslib.org/
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
- RNAstructure: https://rna.urmc.rochester.edu/RNAstructure.html


Install SEISMIC-RNA and its Python dependencies
------------------------------------------------------------------------

After you have installed your chosen non-Python dependencies, run this
command to install SEISMIC-RNA and all its Python dependencies::

    pip install --upgrade seismic-rna

The above command will install the latest version of SEISMIC-RNA and its
dependencies. If you want to install a specific version instead, then
run this command (substituting ``1.2.3`` with the version you want)::

    pip install seismic-rna==1.2.3


.. _bowtie2: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _rnastructure: https://rna.urmc.rochester.edu/RNAstructure.html
.. _samtools: https://www.htslib.org/
