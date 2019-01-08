<p align="center">
  <img src="https://github.com/sivaramambikasaran/HODLR/blob/master/docs/source/images/HODLR.svg" width="250" height="250" alt="Logo of HODLRlib"/>
</p>

# HODLRlib
[![Build Status](https://travis-ci.org/sivaramambikasaran/HODLR.svg?branch=master)](https://travis-ci.org/sivaramambikasaran/HODLR)
[![Documentation Status](https://readthedocs.org/projects/hodlr/badge/?version=latest)](https://hodlr.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/cd427ae7fd414c0cb2a0e1c7d201b2cb)](https://www.codacy.com/app/sivaramambikasaran/HODLR?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=sivaramambikasaran/HODLR&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/sivaramambikasaran/HODLR/badge.svg?branch=master)](https://coveralls.io/github/sivaramambikasaran/HODLR?branch=master)[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Version 3.1415](https://img.shields.io/badge/version-3.1415-brightgreen.svg)](https://github.com/sivaramambikasaran/HODLR)
[![arXiv](https://img.shields.io/badge/math.NA-arXiv%3A1405.0223-%23B31B1B.svg)](https://arxiv.org/abs/1405.0223)
[![star this repo](http://githubbadges.com/star.svg?user=sivaramambikasaran&repo=HODLR&style=flat)](https://github.com/sivaramambikasaran/HODLR)
[![fork this repo](http://githubbadges.com/fork.svg?user=sivaramambikasaran&repo=HODLR&style=flat)](https://github.com/sivaramambikasaran/HODLR/fork)
[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/sivaramambikasaran/HODLR/)
[![PR Welcome](https://img.shields.io/badge/PR-welcome-brightgreen.svg)](http://makeapullrequest.com) 
[![Built by SAFRAN](https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images/built_by_SAFRAN.svg)](http://sivaramambikasaran.com/research/)
[![C++](https://img.shields.io/badge/language-C%2B%2B-brightgreen.svg)](http://www.cplusplus.com/)


HODLRlib is a flexible library for performing matrix operations like matrix-vector products, solving and determinant computation in near-linear complexity(for matrices that resemble a HODLR structure). The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the original articles[[1](https://link.springer.com/article/10.1007/s10915-013-9714-z)][[2](https://arxiv.org/abs/1405.0223)]. Low-rank approximation of the appropriate blocks are obtained using the rook pivoting algorithm. The domain is sub-divided based on a KDTree. The solver is fairly general, works with minimal restrictions and has been parallelized using OpenMP.

For more details on the usage of the library, visit the [documentation](https://hodlr.readthedocs.io/en/latest/) page.

## Features

Fast matrix vector products: Obtains <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//af44b92b9a0ae94e08b5e1e8abce573e.svg?invert_in_darkmode" align=middle width=21.723786149999988pt height=22.465723500000017pt/> at a cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//2d74209c531ea025d06c0a66dbbd0bb1.svg?invert_in_darkmode" align=middle width=85.780695pt height=24.65753399999998pt/>

Fast solution of linear systems: Solves linear systems <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//70681e99f542745bf6a0c56bd4600b39.svg?invert_in_darkmode" align=middle width=50.69621369999999pt height=22.831056599999986pt/> in <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//dad06decfe9b6527d7a6d23885d23d04.svg?invert_in_darkmode" align=middle width=95.43830174999998pt height=29.534320200000014pt/>

Fast symmetric factorization: Compute <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//c0b7ce204101fd85e5ae745c31f7781f.svg?invert_in_darkmode" align=middle width=79.39666349999999pt height=27.6567522pt/> at a cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//dad06decfe9b6527d7a6d23885d23d04.svg?invert_in_darkmode" align=middle width=95.43830174999998pt height=29.534320200000014pt/>

Fast determinant computation: Compute the determinant at an additional cost of <img src="https://cdn.jsdelivr.net/gh/sivaramambikasaran/HODLR@master/docs/source/images//2d74209c531ea025d06c0a66dbbd0bb1.svg?invert_in_darkmode" align=middle width=85.780695pt height=24.65753399999998pt/>

### Version 3.1415

Date: January 6th, 2019

Copyleft 2019: Sivaram Ambikasaran

Developed by Sivaram Ambikasaran, Karan Raj Singh, Shyam Sundar Sankaran

### License

This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/.>
