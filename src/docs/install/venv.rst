
.. _virtual-envs:

Virtual Environments
========================================================================

We *highly* recommend installing SEISMIC-RNA into a virtual environment.
Doing so will likely spare you much frustration well into the future.
If you are not familiar with creating and managing virtual environments,
then we recommend that you read this brief primer.


Why a virtual environment is worth the effort of setting it up
------------------------------------------------------------------------

Virtually every piece of software that you use requires other pieces of
software to be installed. Why? Because software developers want to focus
on implementing their own ideas, rather than "re-inventing the wheel" by
implementing algorithms that have already been written.

Suppose that you are developing a piece of software that requires linear
algebra routines (e.g. matrix multiplication), as does SEISMIC-RNA.
Because these routines are needed ubiquitously, algorithms that perform
them have been distributed in many excellent software packages, such as
`NumPy`_ for Python. Using a pre-made, well-tested package would save
you all the time you would need to plan, code, test, debug, and optimize
your own algorithms (which would take *a lot* of time). Thus, it would
expedite the development of your project to borrow routine tasks from
pre-existing software packages. This borrowing of code is called a
"dependency", and SEISMIC-RNA has dozens.

While dependencies alleviate issues of re-inventing the wheel, they can
also cause problems if dependencies conflict with each other. Suppose
you use two pieces of software regularly at work: call them software A
version 1 (Av1) and software B version 1 (Bv1). Both are dependent on
software C version 1 (Cv1). All goes well until Bv1 becomes deprecated
and you need to update it to Bv2. Bv2 is incompatible with Cv1 and needs
Cv2. But Av1 is incompatible with Cv2, and no upgrade is available that
is compatible with Cv2. So you need to have two versions of software C
installed simultaneously: Cv1 for Av1, and Cv2 for Bv2. How can you
accomplish this feat without mixing up the two versions?


How a virtual environment lets conflicting sofware versions coexist
------------------------------------------------------------------------

Virtual environments let you keep more than one version of a piece of
software on one computer system, without uninstalling and reinstalling
the software every time you needed to use a different version (which
would, of course, be a pain). Instead, you can "activate" a "virtual
environment" that contains the specific version you need. When finished
using that version, you can "deactivate" the virtual environment, which
hides that version of the software.

A virtual environment is simply a collection of software in which each
piece of software has a specific version that you can control, as well
as a list of its dependencies and their versions. The aforementioned
example incompatibility among software A, B, and C could be solved by
installing Av1 and Bv2 into separate virtual environments; call them
`venv1` and `venv2`, respectively. Into `venv1`, you would install Av1
and Cv1; and into `venv2`, you would install Bv2 and Cv2. If you needed
to use softwareA, you would activate `venv1` so that you could run Av1
and Av1 could see Cv1 (but not Cv2). Conversely, when you needed to use
software B, you would activate `venv2` so that you could run Bv2 and Bv2
could see Cv2 (but not Cv1).


.. _venv-environment-managers:

Creating, activating, and deactivating virtual environments
------------------------------------------------------------------------

Creation, deletion, activation, and deactivation of virtual environments
is performed by a piece of software called an "environment manager". Two
popular environment managers for Python are `venv`_ (Python's built-in
virtual environment manager) and `conda`_ (the package manager for the
`Anaconda`_ distribution of Python). If you are unsure of which to use,
then we recommend `conda`_ because it is also a package manager (see
:ref:`venv-package-managers`).

For example, to create a new virtual environment named `seismic` for
SEISMIC-RNA using `conda`, you would enter the following command::

    conda create -n seismic python=3.10

Note that Python 3.10 is currently the only version of Python that is
compatible with SEISMIC-RNA and all of its dependencies.

To activate the environment `seismic`, enter the following command::

    conda activate seismic

Note that it is advisable to give your virtual environments short names
because you will need to type their names every time you activate them.

To deactivate the active virtual environment when you are finished,
enter the following command::

    conda deactivate

Note that you do *not* need to type the name of the environment after
`deactivate` because this command simply deactivates the environment
that is currently active.


.. _venv-package-managers:

Managing software in a virtual environment
------------------------------------------------------------------------

It would be tedious (and risky) to install, uninstall, upgrade, and
downgrade, and all of your software manually. Every time you changed a
piece of software, you would need to install a compatible version of
every one of its dependencies. Fortunately, these tasks can be automated
using another piece of software called a "package manager". (A "package"
is just a specific version of a piece of software and all its code that
you can download.) Two popular package managers for Python are `pip`_
(short for "Package Installer for Python", Python's official package
manager) and `conda`_, which is also an environment manager (see
:ref:`venv-environment-managers`).


.. _Anaconda: https://docs.anaconda.com/free/anaconda/index.html
.. _conda: https://docs.conda.io/en/latest/
.. _NumPy: https://numpy.org/
.. _pip: https://pip.pypa.io/en/stable/
.. _venv: https://docs.python.org/3/library/venv.html
