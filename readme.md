1D HODLR
========================

Version:
========================
1.1.1

Author:
========================
Vaishnavi G.
Contact details: vaishnavihp@gmail.com

Dependencies:
========================
You would need to have g++ and the linear algebra library `Eigen` available @ <http://eigen.tuxfamily.org> installed to run this code.

Files:
========================
HODLR.hpp	:	Header file in C++
Makefile.mk	:	Makefile
readme.md	:	this file
testHODLR.cpp	:	source code in C++

Implementation Details:
========================
This is mainly developed for testing the accuracy and speed of the HODLR. 
This is based on what is discussed in Chapter-7 of the thesis: 'S. Ambikasaran, Fast Algorithms for Dense Numerical Linear Algebra and Applications,
Ph.D. thesis, Stanford University, Palo Alto, CA, 2013'.
An extended sparse matrix is formed which is solved using Eigen's Sparse Solver.
The low rank interaction has been computed using chebyshev interpolation.
Currently, the code is optimized for the kernels 1/(R+0.5), exp(-2.0*R) and log(R+0.5).
You can change the kernels in the makefile by changing KERNEL. 
More kernels can be added by editing the 'getInteraction' function of 'userkernel' class.


Running the code:
========================
To run the code, you would first need to make the Makefile. 

	make -f Makefile.mk

Then we need to run the executable, i.e.,

	./testHODLR L nChebNodes nLevels

where `nLevels` is a positive integer indicating the number of levels in the tree, `nChebNodes` is a positive integer indicating the number of Chebyshev nodes(=rank)  and `L` indicates the semi-length of the simulation rod/line which is from [-L,L]. The number of particles would be `pow(2,nLevels+1)*nChebNodes`.

For instance, below is a sample execution of the code and its output.

It is always good to clean using make clean before running the code, i.e.,
	
	make -f Makefile.mk clean

Then make the file

	make -f Makefile.mk

Run the generated executable as, for instance,

	./testHODLR 1.0 20 8

The generated output will be as follows

	Number of particles is: 10240

	Time taken to create the tree is: 0.00159747

	Time taken to assemble the matrix is: 0.837276

	Time taken to solve is: 3.30924

	Total time taken is: 4.14811

	Apply time taken is: 4.14652

	Total Speed in particles per second is: 2468.59

	Apply Speed in particles per second is: 2469.54

	Performing Error check...

	For Line number: 234

	Error is: 3.8311e-09

	Time taken to compute error is: 0.018475
