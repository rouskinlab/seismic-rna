# A Brief Primer on Virtual Environments


## Why a virtual environment is worth the effort of setting it up

Virtually every piece of software you use depends on other pieces of software.
This strategy enables software developers to avoid "re-inventing the wheel".
Suppose that you are developing a piece of software that requires many linear
algebra routines (e.g. matrix multiplication), as does SEISMIC-RNA. Attempting
to write your own code for these routines would slow down development, bloat the
code base, and very likely produce relatively sluggish software full of bugs.
Much better to borrow linear algebra routines from an existing package that
experts have been developing, debugging, and optimizing for years (e.g. NumPy).
This borrowing of code is termed a dependency. SEISMIC-RNA has dozens of them.

While dependencies alleviate issues of re-inventing the wheel, they can also
cause problems if dependencies conflict with each other. Suppose you use two
pieces of software regularly at work: call them software A version 1 (Av1) and
software B version 1 (Bv1). Both are dependent on software C version 1 (Cv1).
All goes well until Bv1 becomes deprecated and you need to update it to Bv2.
Bv2 is incompatible with Cv1 and requires Cv2. But Av1 is incompatible with Cv2,
and no upgrade is available that is compatible with Cv2. So you need have two
versions of software C installed simultaneously: Cv1 for Av1, and Cv2 for Bv2.
How can you accomplish this feat without mixing up the two versions?


## What a virtual environment is

Virtual environments let you use multiple versions of a piece of software on one
computer system, without needing to uninstall and reinstall the software every
time you needed to use a different version (which would, of course, be a pain).
Instead of tediously installing the version of the software before every use,
you can "activate" a "virtual environment" that contains the specific version
you need. When finished using that version, you can "deactivate" the virtual
environment, thereby hiding that version of the software.

A virtual environment is simply a collection of software, where each piece of
software has a specific version that you can control. The virtual environment
can be "activated" so that you can use the software inside, or "deactivated" so
that the software becomes hidden. Each piece of software also maintains a list
of its dependencies and their versions.

The aforementioned example incompatibility among software A, B, and C could be
solved by installing Av1 and Bv2 into separate virtual environments; call them
`venv1` and `venv2`, respectively. Into `venv1`, you would install Av1 and Cv1;
and into `venv2`, you would install Bv2 and Cv2. When you need to use software
A, you would activate `venv1` so that you could run Av1 and Av1 could see Cv1
(but not Cv2). Conversely, when you need to use software B, you would activate
`venv2` so that you could run Bv2 and Bv2 could see Cv2 (but not Cv1).


## Installing, uninstalling, and versioning software in a virtual environment

It would be extremely tedious (and risky) to install, upgrade, downgrade, and
uninstall your software manually. Every time you changed a piece of software,
you would need to install a compatible version of every one of its dependencies.
Fortunately, there exists software to manage the versions of other software.
Such a piece of software is called a "package manager" (a specific version of a
piece of software, and all its code that you can download, is called a package).

Two very popular package managers for Python are `pip` (an abbreviation of
"Package Installer for Python", which is Python's official package manager:
<https://pip.pypa.io/en/stable/index.html>) and `conda` (the package manager for
the Anaconda distribution of Python: <https://docs.conda.io/en/latest/>).


## Creating, activating, and deactivating virtual environments

While package managers manage the versions of software installed, another kind
of software, an "environment manager", creates, activates, and deactivates the
virtual environments in which software can be installed. The virtual environment
must be created and then activated before any software can be installed into it.

Two very popular environment managers for Python are `venv` (Python's built-in
virtual environment manager, <https://docs.python.org/3/library/venv.html>) and
`conda` (which is both a package manager and an environment manager:
<https://docs.conda.io/en/latest/>).

For example, to create a new virtual environment named `srna` for SEISMIC-RNA
using `conda`, you could enter the following command:

``` conda create -n srna python=3.10 ```

Note that Python 3.10 is currently the only version of Python that is compatible
with SEISMIC-RNA and all of its dependencies.

To activate the environment `srna`, enter the following command:

``` conda activate srna ```

Note that it is advisable to give your virtual environments short names because
you will need to type their names every time you activate them.

To deactivate the active virtual environment when you are finished or want to
switch to a different environment, enter the following command:

``` conda deactivate ```
