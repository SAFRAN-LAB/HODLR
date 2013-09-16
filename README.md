#HODLR_SOLVER  

This is the first public release of the HODLR_SOLVER  library.  
Date: July 31st, 2013

**Version 3.1 - First external release.**

%% Copyleft 2013: Sivaram Ambikasaran
%% Developed by Sivaram Ambikasaran
%% Contact: <siva.1985@gmail.com> (Sivaram)
%%    
%% This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license.      
%% The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not %% distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>  

###FILES:

The following files must be found inside the directory:

1. HODLR_Solver.hpp
2. HODLR_Solver.cpp
3. HODLR_Solver_Input.cpp
4. HODLR_Solver_Tests.cpp
5. makefile_HODLR_Solver.mk
6. makefile_HODLR_Solver_Tests.mk
7. README.md

###SETTING THINGS UP:

1. To run this package, you need to have **Eigen**.

2. Download Eigen from here: <http://eigen.tuxfamily.org/index.php?title=Main_Page>

3. There is a sample input file named "HODLR_Solver_Input.cpp" and another file named "HODLR_Solver_Tests.cpp", which runs the solver for a wide range of system sizes.

4. Go to the directory where Makefile is in, then key in the following command in the terminal:

		make -f makefile_HODLR_Solver.mk

The code should now run.

5. To run some tests for a wide range of system sizes, key in the following command in the terminal:

		make -f makefile_HODLR_Solver_Tests.mk

6. Some default matrices are provided in the HODLR_Solver.hpp. The appropriate matrix can be enabled by changing the makefile

7. Read through the comments in the file HODLR_Solver.hpp, which should be self-explanatory.

8. More details will be added later.
