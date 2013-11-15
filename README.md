#HODLR SOLVER: A blazingly fast direct solver for dense linear systems.

This is an implementation of a blazingly fast direct solver for dense linear systems discussed in the article: "An O(N log (N))  Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices". The solver has been optimized and the running times of the solver and now orders of magnitude faster than the running times reported in the article. The solver has also been extended to matrices not necessarily arising out of kernels. Further, the low-rank approximation of the appropriate blocks are obtained using partial pivoted LU algorithm.

To give a rough idea of the running time, the solver takes 23 seconds (this is time taken from the time you press the enter key on the keyboard to run your code and to get the final result) for a system size of 1 million, where the matrix is of the form

		A	=	s^2 I + B

where

		B(i,j) = sin(R(i,j))/R(i,j)

where R(i,j) is the distance between the points x(i) and x(j), where x(i), x(j) are distributed on a 1D manifold.

<table>
    <tr>
        <td>System size</td> <td>Time taken in seconds</td>
    </tr>
    <tr>
	<td>10 thousand</td> <td>0.2</td>
    </tr>
    <tr>
	<td>100 thousand</td> <td>2.1</td>
    </tr>
    <tr>
	<td>1 Million</td> <td>22.9</td>
    </tr>
</table>

**Author**

Sivaram Ambikasaran <siva.1985@gmail.com>

**Version 3.14**

Date: November 14th, 2013

Copyleft 2013: Sivaram Ambikasaran

Developed by Sivaram Ambikasaran

**License**

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>

###FILES:

The following files must be found inside the directory:

1. HODLR_Node.cpp
2. HODLR_Node.hpp
3. HODLR_Tree.cpp
4. HODLR_Tree.hpp
5. HODLR_Test.cpp
6. get_Matrix.cpp
7. get_Matrix.hpp
8. partial_Piv_LU.cpp
9. partial_Piv_LU.hpp
10. makefile_HODLR_Test.mk
11. README.md

###SETTING THINGS UP:

1. To run this package, you need to have **Eigen**.

2. Download Eigen from here: <http://eigen.tuxfamily.org/index.php?title=Main_Page>

3. There is a sample input file named "HODLR_Test.cpp". This calls the features the code can handle.

4. Go to the directory where Makefile is in, then key in the following command in the terminal:

		make -f makefile_HODLR_Test.mk

5. The makefile provides you with different kernels, you can change. If you want to add more kernels, go and change the files get_Matrix.cpp.

6. Read through the comments in all the files, which should be self-explanatory.

7. More details will be added later.