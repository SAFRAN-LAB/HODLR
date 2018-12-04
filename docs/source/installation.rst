************
Installation
************

Downloading the Source
-----------------------

:math:`\texttt{HODLR}` is distributed using the git version control system, and is hosted on Github. The repository can be cloned using::

    git clone https://github.com/sivaramambikasaran/HODLR.git

Dependencies
-------------

- Eigen Linear Algebra Library(get it `here <https://bitbucket.org/eigen/eigen/>`_)
- (optional) An OpenMP enabled compiler (e.g. gcc4.2 or above) is required to use shared-memory parallelism.
- (optional) MKL libraries(:math:`\texttt{HODLR}` has improved performance when compiled against MKL)

Installation
-------------

First set the environment variable `EIGEN_PATH` to the location of your Eigen installation.
This is needed by the CMake script.

    user@computer HODLR$ export EIGEN_PATH=path/to/eigen/

Optionally: set the enviroment variable `MKLROOT` to take advantage of speedups from MKL.

    user@computer HODLR$ export MKLROOT=path/to/mkl/

Key in the required .cpp to be used as input in CMakeLists.txt. Here you also set the
name of the output executable. Then navigate to your build directory and run `cmake path/to/CMakeLists.txt` and run the generated `Makefile` to get your executable

    user@computer build$ cmake path/to/HODLR/
    user@computer build$ make -j n_threads
    user@computer build$ ./executable
