#HODLR_SOLVER  

Date: November 14th, 2013

**Author**

Sivaram Ambikasaran <siva.1985@gmail.com>

**Version 3.14**

%% Copyleft 2013: Sivaram Ambikasaran

%% Developed by Sivaram Ambikasaran

%% Contact: <siva.1985@gmail.com> (Sivaram)

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

3. There is a sample input file named "HODLR_Test.cpp". This calls the features the code can handle.

4. Go to the directory where Makefile is in, then key in the following command in the terminal:

		make -f makefile_HODLR_Test.mk

5. The makefile provides you with different kernels, you can change. If you want to add more kernels, go and change the files get_Matrix.cpp.

6. Read through the comments in all the files, which should be self-explanatory.

7. More details will be added later.