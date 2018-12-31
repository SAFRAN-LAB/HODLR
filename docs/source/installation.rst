*************************
Installation and Building
*************************

Downloading the Source
-----------------------

:math:`\texttt{HODLRlib}` is distributed using the git version control system, and is hosted on Github. The repository can be cloned using::

    git clone https://github.com/sivaramambikasaran/HODLR.git

Dependencies
-------------

- Eigen Linear Algebra Library(get it `here <https://bitbucket.org/eigen/eigen/>`_)
- (optional) An OpenMP enabled compiler (e.g. gcc4.2 or above) is required to use shared-memory parallelism.
- (optional) MKL libraries(:math:`\texttt{HODLRlib}` has improved performance when compiled against MKL)

Installation
-------------

You can either install :math:`\texttt{HODLRlib}` by using the provided install script provided or manually install and link the needed dependencies.

Install Script
^^^^^^^^^^^^^^

The easiest way to get running is to install the needed dependenceis by running the ``install.sh`` provided in the root level of this repository::

    user@computer HODLR$ source install.sh

The above command should create a folder ``deps/`` in the current directory with the needed dependencies. Additionally, the script should set the environment variables that would be needed during the build and execution stages.

Manually Installing
^^^^^^^^^^^^^^^^^^^

First set the enviroment variable ``HODLR_PATH`` to the root level of this repository. This is needed by some of the routines in plotting of the low rank structure for the specific kernel.

Then, set the environment variable ``EIGEN_PATH`` to the location of your Eigen installation. This is needed by the CMake script.::

    user@computer HODLR$ export EIGEN_PATH=path/to/eigen/

Optionally: set the enviroment variable `MKLROOT` to take advantage of speedups from MKL.::

    user@computer HODLR$ export MKLROOT=path/to/mkl/

Building and Executing
----------------------

Key in the required .cpp to be used as input under ``INPUT_FILE`` in CMakeLists.txt. Here you also set the name of the output executable under ``OUTPUT_EXECUTABLE_NAME``. Then navigate to your build directory and run ``cmake path/to/CMakeLists.txt`` and run the generated ``Makefile`` to get your executable::

    user@computer build$ cmake path/to/HODLR/
    user@computer build$ make -j n_threads
    user@computer build$ ./executable
