=====================
Installation
=====================

.. note::

    If you want to use other modules than ``dreem.draw``, make sure that you have all the dependencies installed. See :ref:`Dependencies` for more information.

Using Docker(Work in progress)
------------
.. note::
    This method of installation assume that docker is installed and running
    

    **There are two options for obtaining and installing the docker image:**

        A. By pulling the image from Docker hub (not recommended at this time)
            ``docker pull rouskin/dreem``
            ``mkdir docker-data``
            
        B. By cloning the repo from github
            ``\n
            git clone https://github.com/rouskinlab/dreem.git\n
            cd dreem\n
            docker build -t dreem . \n
            ``


    **Once done installing and building the image, a container can now be made from that image, in which the analysis can be done:**
        1. First move all the data files that are necessary for the type of analysis that will be done. In most simple cases 
        this will be three files: two fastq files which contains the paired-end reads, and one fasta which contains all the reference sequences.
        
        2. For the docker container to access these, the docker container will mount the dreem/docker-data in the users filesystem, to a folder within the container called var/data/.
        Therefore, for the container to access the desired files, the arguments given to the containers entry point must be given in the form of the path as the container will view it.
        Anything in the docker-data folder will be know to the docker container as var/data/<path/within/docker-data> and therefore the command will be given as:

            ``docker run -iv $(pwd)/docker-data/:/var/data dreem --fasta /var/data/ref.fasta --fastq1 /var/data/sample_R1.fastq --fastq2 /var/data/sample_R2.fastq --temp-dir /var/data/temp/ --out-dir /var/data/out/ --verbose --log /var/data/dreem_docker_run.log``
        
        3. If run correctly, the docker-data folder should be populated with a folder labeled temp, where intermediary files are stored;
        an out folder, where the analysis results are stored; as well as a log file, which details the steps that occured during this executation of the dreem pipeline.

        The containers access to the files can be a little confusing and unintuitive, so below is an example of the correct arguments to be used for the files available in the repo:

            ``docker run -iv $(pwd)/docker-data/:/var/data dreem --fasta /var/data/3509.fasta --fastq1 /var/data/3509-O-flank_1=bi1-ms2-DB_R1.fastq --fastq2 /var/data/3509-O-flank_1=bi1-ms2-DB_R2.fastq --temp-dir /var/data/temp/ --out-dir /var/data/out/ --verbose --log /var/data/dreem_docker_run.log``








Using Pypi and pyenv (Work in progress)
---------------------------------------

::

    [TODO]

.. note::

    [TODO]   

Using Conda
-----------

::

    [TODO]


Using Source
------------------------------------

.. note::
    **In order to install through source code, both conda and x-code (if running on mac OSX).**

    * Instructions for conda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
    
    * Instructions for xcode command line tools: https://www.freecodecamp.org/news/install-xcode-command-line-tools/  
        (this link is not offical and could die at somepoint, but installing xcdoe commandline tools is pretty easy to find with a quick google search)
    
    * Unfortunautly, Windows OS is not friednly to bioinformatics pipelines because of all the 
        dependencies, so to run on windows it is suggested you run on Docker through WSL2. 
        Instructions for installing on WSL2 can be found https://docs.docker.com/desktop/windows/wsl/

::

   cd path/to/where/you/want/dreem
   git clone https://github.com/rouskinlab/dreem.git
   cd dreem
   conda env create -f dreem/env.yml
   conda activate dreem
   pip install .
   pytest 


.. note::


    --the final line ``pytest`` is not required but helpful in ensuring all the correct dependencies have been installed and dreem can access them


