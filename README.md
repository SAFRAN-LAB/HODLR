<p align="center">
  <img src="https://github.com/shyams2/HODLR/blob/master/docs/source/images/HODLR.svg" alt="Logo of HODLRlib"/>
</p>

# HODLRlib
[![Build Status](https://travis-ci.org/sivaramambikasaran/HODLR.svg?branch=master)](https://travis-ci.org/sivaramambikasaran/HODLR)
[![Documentation Status](https://readthedocs.org/projects/hodlr/badge/?version=latest)](https://hodlr.readthedocs.io/en/latest/?badge=latest)


HODLRlib is a flexible library for performing matrix operations like matrix-vector products, solving and determinant computation in near-linear complexity(for matrices that resemble a HODLR structure). The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the original articles[[1](https://link.springer.com/article/10.1007/s10915-013-9714-z)][[2](https://arxiv.org/abs/1405.0223)]. Low-rank approximation of the appropriate blocks are obtained using the rook pivoting algorithm. The domain is sub-divided based on a KDTree. The solver is fairly general, works with minimal restrictions and has been parallelized using OpenMP.

For more details on the usage of the library, visit the [documentation](https://hodlr.readthedocs.io/en/latest/) page.

### Features:

Fast matrix vector products: Obtains <img src="https://cdn.jsdelivr.net/gh/shyams2/HODLR@master/docs/source/images//605193f107171aff38e18841c2fc20a1.svg?invert_in_darkmode" align=middle width=37.24875329999999pt height=22.465723500000017pt/> at a cost of <img src="https://cdn.jsdelivr.net/gh/shyams2/HODLR@master/docs/source/images//2d74209c531ea025d06c0a66dbbd0bb1.svg?invert_in_darkmode" align=middle width=85.780695pt height=24.65753399999998pt/>

Fast solution of linear systems: Solves linear systems <img src="https://cdn.jsdelivr.net/gh/shyams2/HODLR@master/docs/source/images//70681e99f542745bf6a0c56bd4600b39.svg?invert_in_darkmode" align=middle width=50.69621369999999pt height=22.831056599999986pt/> in <img src="https://cdn.jsdelivr.net/gh/shyams2/HODLR@master/docs/source/images//dad06decfe9b6527d7a6d23885d23d04.svg?invert_in_darkmode" align=middle width=95.43830174999998pt height=29.534320200000014pt/>

Fast determinant computation: Compute the determinant at an additional cost of <img src="https://cdn.jsdelivr.net/gh/shyams2/HODLR@master/docs/source/images//2d74209c531ea025d06c0a66dbbd0bb1.svg?invert_in_darkmode" align=middle width=85.780695pt height=24.65753399999998pt/>

#### Version 3.1415

Date: November 29th, 2018

Copyleft 2018: Sivaram Ambikasaran

Developed by Sivaram Ambikasaran, Karan Raj Singh, Shyam Sundar Sankaran

#### License:

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>
