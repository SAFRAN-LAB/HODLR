#HODLR SOLVER: A fast and accurate direct solver for dense linear systems

This is an extension of the fast direct solver discussed in the article: "An O(N log (N)) Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices". The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the article. Low-rank approximation of the appropriate blocks are obtained using partial pivoted LU algorithm. The domain is sub-divided based on a KDTree. The solver is fairly general and works with minimal restrictions.

To give a rough idea of the running time, for a system size of 1 million the solver takes 23 seconds for a 1D problem and 86 seconds for a 2D problem (this is time taken from the time you press the return key on the keyboard to run your code and to get the final result). The matrix is of the form

		A	=	s^2 I + B

where B(i,j) depends on R(i,j), the distance between the points x(i) and x(j). The computed answer is accurate upto 12-13 digits. It is to be noted that even for higher dimensions (2D and 3D), the solver is still way faster than the conventional direct solvers, though the scaling may no longer be almost linear. The scaling depends on the smoothness of the kernel near the diagonal of the matrix. For instance, if the kernel is log(R) the scaling is almost linear in 1D as shown by the timing below.

***In 1D***
<table>
    <tr>
        <td>System size</td> <td>Time taken</td> <td>Accuracy</td>
    </tr>
    <tr>
	<td>10 thousand</td> <td>0.5 seconds</td> <td>12 digits</td>
    </tr>
    <tr>
	<td>100 thousand</td> <td>6.3 seconds</td> <td>11 digits</td>
    </tr>
    <tr>
	<td>500 thousand</td> <td>40 seconds</td> <td>11 digits</td>
    </tr>
</table>

***In 2D***
<table>
    <tr>
        <td>System size</td> <td>Time taken</td> <td>Accuracy</td>
    </tr>
    <tr>
	<td>10 thousand</td> <td>6.8 seconds</td> <td>12 digits</td>
    </tr>
    <tr>
	<td>100 thousand</td> <td>390 seconds</td> <td>11 digits</td>
    </tr>
</table>


If the kernel is sin(R)/R (smooth close to the diagonal), the scaling is linear in all dimensions.

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
    <tr>
	<td>500 thousand</td> <td>115 seconds</td> <td>12 digits</td>
    </tr>
</table>

**Author**

Sivaram Ambikasaran <siva.1985@gmail.com>

**Citation**

If you use the implementation or any part of the implementation in your work, kindly cite as follows:

***Article***

@article{ambikasaran2013fastdirect,</br>
  author={{A}mbikasaran, {S}ivaram and {D}arve, {E}ric},</br>
  title={An $\mathcal{O}(N \log N)$ Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices},</br>
  journal={Journal of Scientific Computing},</br>
  year={2013},</br>
  volume={57},</br>
  number={3},</br>
  pages={477--501},</br>
  month={December},</br>
  publisher={Springer}
}

***Code***

@MISC{ambikasaran2013HODLR,</br>
  author = {{A}mbikasaran, {S}ivaram},</br>
  title = {A fast direct solver for dense linear systems},</br>
  howpublished = {https://github.com/sivaramambikasaran/HODLR_Solver},</br>
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

###SETTING THINGS UP:

1. To run this package, you need to have **Eigen**.

2. Download Eigen from here: <http://eigen.tuxfamily.org/index.php?title=Main_Page>

3. There is a sample input file named "HODLR_Test.cpp" in the directory './examples/'. This calls the features the code can handle.

4. Go to the directory where makefile is in, then key in the following command in the terminal:

		make -f makefile.mk

5. Once your run the make command, the executables are created in the directory named './exec/'. To run the code, go into the 'exec' directory and call './HODLR_Test'.

5. You can change the kernels in the makefile by changing KERNEL. More kernels can be added by editing the function

		double get_Matrix_Entry(const unsigned i, const unsigned j)

 in the file

		get_Matrix.cpp

6. The dimension of can be changed by changing DIM in the makefile.

7. Read through the comments in all the files. Most of the function/class/method/variable names are self-explanatory. Read the file HODLR_Test.cpp to understand how to assemble, factor, solver a new HODLR system.

8. More details will be added later.