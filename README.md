#HODLR: Fast direct solver and determinant computation for dense linear systems

This is an extension of the fast direct solver discussed in the article: "An O(N log (N)) Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices". The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the article. Low-rank approximation of the appropriate blocks are obtained using partial pivoted LU algorithm. The domain is sub-divided based on a KDTree. The solver is fairly general and works with minimal restrictions.

To give a rough idea of the running time, for a system size of 1 million the solver takes 23 seconds for a 1D problem and 86 seconds for a 2D problem (this is time taken from the time you press the return key on the keyboard to run your code and to get the final result). The matrix is of the form

		A	=	s^2 I + B

where B(i,j) depends on R(i,j), the distance between the points x(i) and x(j). The computed answer is accurate upto more than 10 digits. It is to be noted that even for higher dimensions (2D and 3D), the solver is still way faster than the conventional direct solvers, though the scaling may no longer be almost linear. The scaling depends on the smoothness of the kernel near the diagonal of the matrix.

**Author**

Sivaram Ambikasaran <siva.1985@gmail.com>

**Citation**

If you use the implementation or any part of the implementation in your work, kindly cite as follows:

***Article***

@article{ambikasaran2013fast,

  author={{A}mbikasaran, {S}ivaram and {D}arve, {E}ric},
  
  title={An $\mathcal{O}(N \log N)$ Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices},
  
  journal={Journal of Scientific Computing},
  
  year={2013},
  
  volume={57},
  
  number={3},
  
  pages={477--501},
  
  month={December},
  
  publisher={Springer}
  
}

***Code***

@MISC{ambikasaran2013HODLR,

  author = {{A}mbikasaran, {S}ivaram},
  
  title = {A fast direct solver for dense linear systems},
  
  howpublished = {https://github.com/sivaramambikasaran/HODLR},
  
  year = {2013}
  
 }

**Version 3.14**

Date: November 14th, 2013

Copyleft 2013: Sivaram Ambikasaran

Developed by Sivaram Ambikasaran

**License**

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>

###DIRECTORIES AND FILES


	./examples/		:	Example input C++ codes; Needed to read input from user or from input file.
	./src/			:	Source code in C++
	./header/		:	Relevant header files
	./exec/			:	Executables for HODLR
	./README.md		:	This file
	./LICENSE.md	:	License file
	./makefile.mk	:	Makefile

##Usage

###DEPENDENCIES:

To run this package, you need to have **Eigen**. If you don't already have it,
download and install Eigen following the instructions
[here](http://eigen.tuxfamily.org/index.php?title=Main_Page).

###BUILD USING CMAKE:

The easiest way to build this library is using [CMake](http://cmake.org/).
In the project directory, run:
```
mkdir build
cd build
cmake ..
make
make test
[sudo] make install # optional
```
this will build the static `hodlr` library and run a few tests. If your
version of the Eigen headers is installed in a non-standard place, you can
change the `cmake` line to:
```
cmake .. -DEIGEN_INCLUDE_DIR_HINTS=/path/to/eigen
```

Your code should include `get_Matrix.hpp` and implement the function
```
double get_Matrix_Entry (const unsigned i, const unsigned j)
```


###CUSTOM MAKEFILE:

1. There is a sample input file named "HODLR_Test.cpp" in the directory './examples/'. This calls the features the code can handle.

2. Go to the directory where makefile is in, then key in the following command in the terminal:

		make -f makefile.mk

3. Once your run the make command, the executables are created in the directory named './exec/'. To run the code, go into the 'exec' directory and call './HODLR_Test'.

4. You can change the kernels in the makefile by changing KERNEL. More kernels can be added by editing the function

		double get_Matrix_Entry(const unsigned i, const unsigned j)

 in the file

		get_Matrix.cpp

5. The dimension of the problem can be changed by changing DIM in the makefile.

6. Read through the comments in all the files. Most of the function/class/method/variable names are self-explanatory. Read the file HODLR_Test.cpp to understand how to assemble, factor, solver a new HODLR system.
