
.. _install-prereqs:

Prerequisites: Python, pip, and Conda
========================================================================

This section describes how to install Python, pip, and Conda if
if they are not already installed on your system (or if you are not sure
whether they are installed). SEISMIC-RNA does not require that pip
and Conda are installed, but we recommend using them because they
make installing SEISMIC-RNA and its dependencies much easier, especially
for novice users.


Install Conda
------------------------------------------------------------------------

Installing Conda installs Python and pip automatically, so it is
the easiest option to start. Visit the webpage for installing Conda
(https://docs.conda.io/projects/conda/en/latest/user-guide/install/),
click the link for your operating system under "Regular installation",
download the installer named "Miniconda installer for ...", and follow
the rest of the instructions. You may answer "yes" to all the default
options if you so desire. **Make sure that you remember the path where
conda is being installed**. On Linux and macOS, the default is
``~/miniconda3``.


Initialize Conda
------------------------------------------------------------------------

After installing Conda, you will need to initialize it so that it is
set up to run properly. Initialize by typing this in the command line,
making sure to replace ``/path/to/your/conda`` with the actual path in
which Conda was installed, plus ``/condabin/conda`` (for example, if
the installation path were ``~/miniconda3``, then the path to type would
be ``~/miniconda3/condabin/conda``)::

    /path/to/your/conda init

If it worked, then you will see a message that looks similar to this::

    no change     /home/myname/conda/condabin/conda
    no change     /home/myname/conda/bin/conda
    no change     /home/myname/conda/bin/conda-env
    no change     /home/myname/conda/bin/activate
    no change     /home/myname/conda/bin/deactivate
    no change     /home/myname/conda/etc/profile.d/conda.sh
    no change     /home/myname/conda/etc/fish/conf.d/conda.fish
    no change     /home/myname/conda/shell/condabin/Conda.psm1
    no change     /home/myname/conda/shell/condabin/conda-hook.ps1
    no change     /home/myname/conda/lib/python3.9/site-packages/xontrib/conda.xsh
    no change     /home/myname/conda/etc/profile.d/conda.csh
    changed       /home/myname/.bash_profile

Then, **close and reopen** your command line (the changes will not take
effect until you restart the command line).


Check your installation
------------------------------------------------------------------------

Verify that Conda is initialized by typing the following **after**
you have restarted the command line::

    which conda

If Conda is initialized, then this command will print either the
path to the Conda executable, such as this::

    /home/myname/miniconda3/condabin/conda

or a several lines of code that look like this::

    conda () {
        \local cmd="${1-__missing__}"
        case "$cmd" in
            (activate | deactivate) __conda_activate "$@" ;;
            (install | update | upgrade | remove | uninstall) __conda_exe "$@" || \return
                __conda_reactivate ;;
            (*) __conda_exe "$@" ;;
        esac
    }

If it is not initialized, then you will get an error message that says
Conda was not found.


Troubleshooting Conda installation/initialization
------------------------------------------------------------------------

If you encounter any errors, refer to the Conda documentation at
https://docs.conda.io/projects/conda/en/latest/user-guide/troubleshooting.html.
