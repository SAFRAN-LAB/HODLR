##FORTRAN Wrappers for the HODLR Library

This folder comprises a set of C and FORTRAN files
which demonstrate a FORTRAN driver calling several of 
the core routines from the HODLR library.

**Author**

* FORTRAN Wrappers Travis Askham <askhamwhat@gmail.com>

* HODLR Library Sivaram Ambikasaran

**Citation**

See README.md of parent HODLR/ directory.

**License**

The files contained in THIS folder (the FORTRAN wrappers) are 
available under the MIT License.

The MIT License (MIT)

Copyright (c) 2014 Travis Askham

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


###DIRECTORIES AND FILES

* HODLR_Extern.cpp - a C++ class to aid interoperability with FORTRAN
* Matrix_Entry_Routine.cpp - contains the function used to evaluate the
matrix entries (should be edited based on application)
* driver2.f - a FORTRAN driver file. Contains the interfaces for the 
wrapper routines and demonstrates usage in the context of RBF 
interpolation in 1D.
* fort_hodlr_wrappers.cpp - a set of C++ routines which are callable 
from FORTRAN
* fort_hodlr_wrappers.hpp - extern "C" wrapper for the C++ routines,
prevents the linker from mangling the name of the subroutines.
* driver2.mk - makefile. See below.

###CAPABILITIES

With these files, one can (from FORTRAN):

* Initialize an instance of the HODLR_Tree class (needed to perform
the HODLR calculations), corresponding to the system matrix.
* Factor this matrix.
* Solve a system Ax = b.
* Perform a matrix multiplication b = Ax.
* Compute the log-determinant of this matrix.

See fort_hodlr_wrappers.cpp and driver2.f for usage. The routines
have a reasonable amount of documentation.

###COMPILE

* In this folder, enter

> make -f driver2.mk

to compile. You may have to change the location of Eigen/Dense 
in the makefile. You can change the compilers as well, though it
has only been tested with gfortran and g++.

* To run

> ./DRIVER2

###DEPENDENCIES:

These routines rely on the HODLR and Eigen libraries.

###SOME NOTES ON USAGE

The system considered is of the form

> A = diag + B.

The envisioned usage of the routines is as follows:

1. The FORTRAN driver determines a set of data which defines
the entries of the system matrix. 
2. This data is then stored in the diag, xyz, xyztarg, and work 
arrays in any manner that makes sense to the user.
3. When the instance of the HODLR_Tree is created (the system
matrix) this data is passed to the routine which calculates
the matrix entries.
4. To change the way the matrix entries are calculated using
these arrays, one can edit the Matrix_Entry_Routine.cpp file.
5. Once the matrix is created, it should be factored using the
supplied routine.
6. After factorization, any of the multiplication, solve, 
and determinant routines may be called.

Note: the system matrix is stored as a pointer to a pointer
to the HODLR_Tree object that stores it. On the FORTRAN side,
this pointer should only be manipulated by using the 
fort_hodlr_wrappers.cpp routines.

When writing your own FORTRAN driver, copy the "interface"
statements and the iso_c_bindings statement from the supplied 
driver (driver2.f) to your file (somewhere between "program" and
data declarations)