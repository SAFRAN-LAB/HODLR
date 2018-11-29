# HODLR: Fast direct solver and determinant computation for dense linear systems
[![Build Status](https://travis-ci.com/shyams2/HODLR.svg?branch=master)](https://travis-ci.com/shyams2/HODLR)

This is an extension of the fast direct solver discussed in the article: "An O(N log (N)) Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices". The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the article. Low-rank approximation of the appropriate blocks are obtained using partial pivoted LU algorithm. The domain is sub-divided based on a KDTree. The solver is fairly general, works with minimal restrictions and has been parallelized using OpenMP. The key features of the solver include

### Features:

```
Fast matrix vector products: Obtains A*x at a cost of O(N log N)
Fast solution of linear systems: Solves linear systems Ax = b in O(N log^2N)
Fast determinant computation: Compute the determinant at an additional cost of O(N log N)
```

### Dependencies:

To run this package, you need to have **Eigen**. If you don't already have it, download and install Eigen following the instructions [here](http://eigen.tuxfamily.org/index.php?title=Main_Page).

### Getting Started:

1. First set the environment variable `EIGEN_PATH` to the location of your Eigen installation. This is needed by the CMake script.

2. Key in the required `.cpp` to be used as input in `CMakeLists.txt`. Here you also set the name of the output executable.

3. There is a sample example file named: "testHODLR.cpp" in the directory './examples/'. This calls the features of the code and reports the timings and the errors. For this demonstration, we assume that `testHODLR.cpp` is given as the input in `CMakeLists.txt` with the output executable named as `example`.

4. For the sake of demonstration, we will be creating the build folder in the current directory itself. This is done using:

    ```
    mkdir build
    cd build
    ```

5. Now in this build folder generate the appropriate `Makefile` by running `cmake ..`. You can now build the system using this `Makefile` to produce the executable:

    ```
    make -j n_threads
    ```

6. Once your run the make command, the executable `example` should be created in the current directory. To run the code, key in:

    ```
    ./example N M d tol
    ```

where `N` is the size of the system you like to handle, `M` is the size of the smallest system you can handle without the fast code (essentially 'M' is the size of the matrix at the leaf nodes), `d` is the dimensionality of the problem considered.`tol` is used to set the tolerance of the computation as `10^{-tol}`
        
Read the file `testHODLR.cpp` to understand how to assemble, factor, solve, and evaluate the log-determinant of a HODLR system. In this example, the matrix assembled is using the gaussian kernel. You can introduce your own kernel function by changing the `getMatrixEntry` function.

#### Version 3.1415

Date: November 29th, 2018

Copyleft 2018: Sivaram Ambikasaran

Developed by Sivaram Ambikasaran, Karan Raj Singh, Shyam Sundar Sankaran

#### License:

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>
