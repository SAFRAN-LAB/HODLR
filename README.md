# HODLR: Fast direct solver, symmetric factorization and determinant computation for dense linear systems

This is an extension of the fast direct solver discussed in the article: "An O(N log (N)) Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices". The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the article. Low-rank approximation of the appropriate blocks are obtained using partial pivoted LU algorithm. The domain is sub-divided based on a KDTree. The solver is fairly general, works with minimal restrictions and has been parallelized using openmp. The key features of the solver include

### FEATURES
	
	Fast matrix vector products: Obtains A*x at a cost of O(N log N)

	Fast solution of linear systems: Solves linear systems Ax=b in O(N log^2N)

	Fast determinant computation: Compute the determinant at an additional cost of O(N log N)

	Fast symmetric factorization: Obtains the symmetric factorization of a symmetric positive definite matrix, i.e., A = WW^T, at a cost of O(N log^2 N)

	Fast symmetric factorization: Solves the symmetric positive definite matrix, at a cost of O(N log N)

	Fast application of symmetric factor to a vector: Obtains W*x at a cost of O(N log N)


### DIRECTORIES AND FILES


	./examples/		:	Example input C++ codes; Needed to read input from user or from input file.
	./src/			:	Source code in C++
	./header/		:	Relevant header files
	./exec/			:	Executables for HODLR
	./README.md		:	This file
	./LICENSE.md	:	License file
	./makefile.mk	:	Makefile

### DEPENDENCIES:

To run this package, you need to have **Eigen**. If you don't already have it, download and install Eigen following the instructions [here](http://eigen.tuxfamily.org/index.php?title=Main_Page).

### CUSTOM BUILD:

1. There is a sample input file named: "testHODLR.cpp" in the directory './examples/'.

2. This calls the features of the code and reports the timings and the errors.

3. Go to the directory where the makefile is and key in the following command in the terminal:
	make -f makefile.mk

4. Once your run the make command, the executables are created in the directory named './exec/'. To run the code, key in

		./exec/HODLR_Test N M d

where 'N' is the size of the system you like to handle, 'M' is the size of the smallest system you can handle without the fast code (essentially M is the size of the matrix at the leaf nodes) and 'd' is the (rough) number of digits of accuracy you need. For instance,
		
		./exec/HODLR_Test 100000 200 12

More kernels can be added by editing the function

		double get_Matrix_Entry(const unsigned i, const unsigned j)

 in the file

		HODLR_Matrix.hpp

Read through the comments in all the files. Most of the function/class/method/variable names are self-explanatory. Read the file HODLR_Test.cpp to understand how to assemble, factor, solve, symmetric factorization, determinant computation of a new HODLR system.

#### Version3.141

Date: November 9th, 2016

Copyleft 2016: Sivaram Ambikasaran

Developed by Sivaram Ambikasaran, Karan Raj Singh

#### License
This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>