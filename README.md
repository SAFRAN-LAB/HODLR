#HODLR SOLVER: Blazingly fast and accurate direct solver for dense linear systems

This is an implementation of a blazingly fast direct solver for dense linear systems discussed in the article: "An O(N log (N)) Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices". The solver has been optimized and the running times of the solver is now massively (a few orders of magnitude) faster than the running times reported in the article. The solver has also been extended to matrices not necessarily arising out of kernels. Low-rank approximation of the appropriate blocks are obtained using partial pivoted LU algorithm. The solver is highly general and works for systems, where the off-diagonal blocks can be efficiently represented as a low-rank matrix in a hierarchical fashion.

To give a rough idea of the running time, the solver takes 23 seconds (this is time taken from the time you press the return key on the keyboard to run your code and to get the final result) for a system size of 1 million. The matrix is of the form

		A	=	s^2 I + B

where

		B(i,j) = sin(R(i,j))/R(i,j)

where R(i,j) is the distance between the points x(i) and x(j). The computed answer is accurate upto 13 digits. It is to be noted that even for higher dimensions (2D and 3D), the solver is still way faster than the conventional direct solvers, though the scaling may no longer be almost linear. The scaling depends on the smoothness of the kernel near the diagonal of the matrix. For instance, if the kernel is sin(R)/R, the timings indicate that the scaling is linear in all dimensions.

***In 1D***
<table>
    <tr>
        <td>System size</td> <td>Time taken</td> <td>Accuracy</td>
    </tr>
    <tr>
	<td>10 thousand</td> <td>0.2 seconds</td> <td>13 digits</td>
    </tr>
    <tr>
	<td>100 thousand</td> <td>2.1 seconds</td> <td>13 digits</td>
    </tr>
    <tr>
	<td>1 Million</td> <td>22.9 seconds</td> <td>13 digits</td>
    </tr>
</table>

***In 2D***
<table>
    <tr>
        <td>System size</td> <td>Time taken</td> <td>Accuracy</td>
    </tr>
    <tr>
	<td>10 thousand</td> <td>0.4 seconds</td> <td>13 digits</td>
    </tr>
    <tr>
	<td>100 thousand</td> <td>5 seconds</td> <td>13 digits</td>
    </tr>
    <tr>
	<td>500 thousand</td> <td>31 seconds</td> <td>12 digits</td>
    </tr>
    <tr>
	<td>1 Million</td> <td>86 seconds</td> <td>12 digits</td>
    </tr>
</table>

***In 3D***
<table>
    <tr>
        <td>System size</td> <td>Time taken</td> <td>Accuracy</td>
    </tr>
    <tr>
	<td>10 thousand</td> <td>1.2 seconds</td> <td>13 digits</td>
    </tr>
    <tr>
	<td>50 thousand</td> <td>6 seconds</td> <td>13 digits</td>
    </tr>
    <tr>
	<td>100 thousand</td> <td>12.5 seconds</td> <td>12 digits</td>
    </tr>
    <tr>
	<td>200 thousand</td> <td>26 seconds</td> <td>12 digits</td>
    </tr>
    <tr>
	<td>400 thousand</td> <td>72 seconds</td> <td>12 digits</td>
    </tr>
</table>

**Author**

Sivaram Ambikasaran <siva.1985@gmail.com>

**Citation**

If you use this as part of your work, kindly cite as follows:

@article{ambikasaran2013fastdirect,
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

@MISC{ambikasaran2013HODLR,
  author = {{A}mbikasaran, {S}ivaram},
  title = {Blazingly fast direct solver for dense linear systems},
  howpublished = {https://github.com/sivaramambikasaran/HODLR_Solver},
  year = {2013}
 }

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

4. Go to the directory where makefile is in, then key in the following command in the terminal:

		make -f makefile_HODLR_Test.mk

5. You can change the kernels in the makefile. More kernels can be added by editing the function

		double get_Matrix_Entry(const unsigned i, const unsigned j)

 in the file

		get_Matrix.cpp

6. Read through the comments in all the files. Most of the function/class/method/variable names are self-explanatory.

7. More details will be added later.