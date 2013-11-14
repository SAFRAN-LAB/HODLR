#HODLR_SOLVER  

Date: November 14th, 2013

**Version 3.14**

%% Copyleft 2013: Sivaram Ambikasaran
%% Developed by Sivaram Ambikasaran
%% Contact: <siva.1985@gmail.com> (Sivaram)
%%    
%% This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license.      
%% The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not %% distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>  

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

3. There is a sample input file named "HODLR_Solver_Input.cpp" and another file named "HODLR_Solver_Tests.cpp", which runs the solver for a wide range of system sizes.

4. Go to the directory where Makefile is in, then key in the following command in the terminal:

		make -f makefile_HODLR_Test.mk

5. To change the matrices, go and change the files get_Matrix.cpp.

7. Read through the comments in all the files, which should be self-explanatory.

8. More details will be added later.