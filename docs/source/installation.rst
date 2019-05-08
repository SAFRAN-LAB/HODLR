*************************
Installation and Building
*************************

Downloading the Source
-----------------------

:math:`\texttt{HODLRlib}` is distributed using the git version control system, and is hosted on Github. The repository can be cloned using::

    git clone https://github.com/sivaramambikasaran/HODLR.git --recursive

The ``--recursive`` flag is argument to ensure that all submodules are also checked out. The Python wrappers for the library make use of ``pybind11`` and would need to be checked out for Python functionality.

Dependencies
-------------

- Eigen Linear Algebra Library (get it `here <https://bitbucket.org/eigen/eigen/>`_)
- (optional) An OpenMP enabled compiler (e.g. gcc4.2 or above) is required to use shared-memory parallelism.
- (optional) MKL libraries (:math:`\texttt{HODLRlib}` has improved performance when compiled against MKL)

**NOTE**: On MacOS, the default compiler is `clang` which doesn't have OpenMP support. You will have to use g++ to make use of the speedups from OpenMP::

    user@computer HODLR$ brew install g++-8
    user@computer HODLR$ export CXX=g++

Installation
-------------

You can either install :math:`\texttt{HODLRlib}` by using the provided install script provided or manually install and link the needed dependencies.

Install Script
^^^^^^^^^^^^^^

The easiest way to get running is to install the needed dependencies by running the ``install.sh`` provided in the root level of this repository::

    user@computer HODLR$ ./install.sh

The above command should create a folder ``deps/`` in the current directory with the needed dependencies. Additionally, the script should set the environment variables that would be needed during the build and execution stages. This only needs to be done once since the environment variables are automatically written to ``.bash_profile``.

Manually Installing
^^^^^^^^^^^^^^^^^^^

First set the environment variable ``HODLR_PATH`` to the root level of this repository. This is needed by some of the routines in the plotting of the low-rank structure for the specific kernel. (**NOTE**: The plotting is carried out using python, and requires the  `matplotlib <https://matplotlib.org/>`_ package to be installed in your python environment)

Then, set the environment variable ``EIGEN_PATH`` to the location of your Eigen installation. This is needed by the CMake script.::

    user@computer HODLR$ export EIGEN_PATH=path/to/eigen/

Optionally: set the environment variable ``MKLROOT`` to take advantage of speedups from MKL.::

    user@computer HODLR$ export MKLROOT=path/to/mkl/

Python Installation
^^^^^^^^^^^^^^^^^^^

Get the dependencies and set environment variables needed for the library through the Install Script of Manual Install routes. Once these have been setup, navigate to the ``python`` subdirectory from the root level of the repository. The Python wrapper for the library can be installed by running::
    
    user@computer python$ python setup.py install 

Testing
-------

C++
^^^

Now, we need to ensure that all the functions of the libraries function as intended. For this purpose, we will be running the script ``test/test_HODLR.cpp``. By default, during a build this file under ``test/`` gets compiled, and would show up under the ``test/`` directory in your build folder. To check this on your computer, run the following lines::

    user@computer HODLR$ mkdir build && cd build
    user@computer build$ cmake ..
    user@computer build$ ./test/test_HODLR

For a succesful test, the final line of output for this run would read:"Reached End of Test File Successfully! All functions work as intended!".

Fortran
^^^^^^^

In order to ensure that the Fortran wrapper works as intended, run the Fortran file found under ``fortran/test/test_HODLR.f90``. In order to run this and ensure correct functioning, type the following at the commandline::

    user@computer HODLR$ cd fortran
    user@computer fortran$ mkdir build && cd build
    user@computer build$ cmake .. -DCMAKE_BUILD_TYPE=TESTING
    user@computer build$ ./test_HODLR

For a succesful test, the final line of output for this run would read:"Reached End of Test File Successfully! All functions work as intended!".

Python
^^^^^^

In order to ensure that the Python wrapper works as intended, run the Python file found under ``python/test/test_HODLR.py``. This can be executed by running the following commands from the terminal::

    user@computer HODLR$ cd python/test
    user@computer test$ python test_HODLR.py

For a succesful test, the final line of output for this run would read:"Reached End of Test File Successfully! All functions work as intended!".

Building and Executing
----------------------

C++
^^^

Key in the required ``.cpp`` to be used as input under ``INPUT_FILE`` in ``HODLRlib/CMakeLists.txt``. Here you also set the name of the output executable under ``OUTPUT_EXECUTABLE_NAME``. Then navigate to your build directory and run ``cmake path/to/CMakeLists.txt`` and run the generated ``Makefile`` to get your executable::

    user@computer build$ cmake path/to/HODLR/
    user@computer build$ make -j n_threads
    user@computer build$ ./executable

Fortran
^^^^^^^

The compilation process for the fortran wrappers is similar to the process for C++: Key in the required ``.f90`` to be used as input under ``INPUT_FILE`` in ``HODLRlib/fortran/CMakeLists.txt``. Here you also set the name of the output executable under ``OUTPUT_EXECUTABLE_NAME``. Then navigate to your build directory and run ``cmake path/to/CMakeLists.txt`` and run the generated ``Makefile`` to get your executable::

    user@computer build$ cmake path/to/HODLR/fortran/
    user@computer build$ make -j n_threads
    user@computer build$ ./executable

Python
^^^^^^

The process of running Python files that make use of ``pyhodlrlib`` follow the standard procedure used to run a python file::

    user@computer example$ python example.py
