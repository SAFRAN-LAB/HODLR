<p align="center">
  <img src="https://github.com/sivaramambikasaran/HODLR/blob/master/docs/source/images/HODLR.svg" width="250" height="250" alt="Logo of HODLRlib"/>
</p>

# HODLRlib

[![Documentation Status](https://readthedocs.org/projects/hodlrlib/badge/?version=latest)](https://hodlrlib.readthedocs.io/en/latest/?badge=latest)

[![C++](https://img.shields.io/badge/language-C%2B%2B-brightgreen.svg)](http://www.cplusplus.com/)
[![Build Status](https://travis-ci.org/sivaramambikasaran/HODLR.svg?branch=master)](https://travis-ci.org/sivaramambikasaran/HODLR)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/cd427ae7fd414c0cb2a0e1c7d201b2cb)](https://www.codacy.com/app/sivaramambikasaran/HODLR?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=sivaramambikasaran/HODLR&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/sivaramambikasaran/HODLR/badge.svg?branch=master)](https://coveralls.io/github/sivaramambikasaran/HODLR?branch=master)

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)

[![Version 3.1415](https://img.shields.io/badge/version-3.1415-brightgreen.svg)](https://github.com/sivaramambikasaran/HODLR)
[![arXiv](https://img.shields.io/badge/math.NA-arXiv%3A1405.0223-%23B31B1B.svg)](https://arxiv.org/abs/1405.0223)

[![star this repo](http://githubbadges.com/star.svg?user=sivaramambikasaran&repo=HODLR&style=flat)](https://github.com/sivaramambikasaran/HODLR)
[![fork this repo](http://githubbadges.com/fork.svg?user=sivaramambikasaran&repo=HODLR&style=flat)](https://github.com/sivaramambikasaran/HODLR/fork)

[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/sivaramambikasaran/HODLR/)
[![PR Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](http://makeapullrequest.com) 

[![Built by SAFRAN](https://img.shields.io/badge/built%20by-SAFRAN-orange.svg)](http://sivaramambikasaran.com/research/)

[![DOI](https://zenodo.org/badge/12858603.svg)](https://zenodo.org/badge/latestdoi/12858603)

[![DOI](http://joss.theoj.org/papers/10.21105/joss.01167/status.svg)](https://doi.org/10.21105/joss.01167)

HODLRlib is a flexible library for performing matrix operations like matrix-vector products, solving and determinant computation in near-linear complexity(for matrices that resemble a HODLR structure). The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the original articles[[1](https://link.springer.com/article/10.1007/s10915-013-9714-z)][[2](https://arxiv.org/abs/1405.0223)]. Low-rank approximation of the appropriate blocks are obtained using the rook pivoting algorithm. The domain is sub-divided based on a KDTree. The solver is fairly general, works with minimal restrictions and has been parallelized using OpenMP.

For more details on the usage of the library, visit the [documentation](https://hodlrlib.readthedocs.io/) page.

## Features

MatVecs: Obtains <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//af44b92b9a0ae94e08b5e1e8abce573e.svg?invert_in_darkmode" align=middle width=21.723786149999988pt height=22.465723500000017pt/> at a cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//a905df5a2fee5cc61be08bca001d96bc.svg?invert_in_darkmode" align=middle width=85.780695pt height=24.65753399999998pt/>

Factorization: Factors the matrix <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode" align=middle width=12.32879834999999pt height=22.465723500000017pt/> into the desired form at a cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//b968ed4db2c93b9d0e799f4fc7300fed.svg?invert_in_darkmode" align=middle width=108.22370954999998pt height=29.534320200000014pt/>

Cholesky-like Symmetric Factorization: Obtains <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//9fe6a39ffcbf0fd1960b9767d054cf6e.svg?invert_in_darkmode" align=middle width=79.39666349999999pt height=27.6567522pt/> at a cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//b968ed4db2c93b9d0e799f4fc7300fed.svg?invert_in_darkmode" align=middle width=108.22370954999998pt height=29.534320200000014pt/>

Solve: Solves linear systems <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//66a8a0c17c80a313cb880fcc6d6392f3.svg?invert_in_darkmode" align=middle width=50.69621369999999pt height=22.831056599999986pt/> at an additional cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//053e1c7d38b655ef637f98b669d34798.svg?invert_in_darkmode" align=middle width=98.5661028pt height=24.65753399999998pt/>

Determinant Computation: Additional Cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//87a88fd17efcaab84de7605c60cd4528.svg?invert_in_darkmode" align=middle width=85.780695pt height=24.65753399999998pt/>

### Version 3.1415

Date: January 6th, 2019

Copyleft 2019: Sivaram Ambikasaran

Developed by Sivaram Ambikasaran, Karan Raj Singh, Shyam Sundar Sankaran

### License

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>
