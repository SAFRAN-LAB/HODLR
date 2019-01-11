---
title: 'HODLRlib: A Library for Hierarchical Matrices'
tags:
  - Hierachical Matrix
  - Fast Algorithms
  - Covariance Matrices
  - Gaussian Processes
  - Interpolation
  - Brownian Dynamics
authors:
  - name: Sivaram Ambikasaran
    affiliation: 1
  - name: Karan Raj Singh
    affiliation: 2
  - name: Shyam Sundar Sankaran
    affiliation: 1
affiliations:
  - name: Department of Mathematics, Indian Institute of Technology Madras
    index: 1
  - name: Department of Computational & Data Sciences, Indian Institute of Science
    index: 2
date: 11th January 2019
bibliography: paper.bib
---

# Summary

``HODLRlib`` is a flexible library for working with matrices that have a Hierarchical Off-Diagonal Low-Rank (HODLR) structure. The current version performs matrix operations like matrix-vector products, solving linear systems, Cholesky-like symmetric factorization and determinant computation in almost linear complexity. A key motivation for ``HODLRlib`` is to provide a high-performance, easy to use library for working with matrices that possess a HODLR structure. 

``HODLRlib`` is an optimized implementation of the ideas illustrated in these articles[@ambikasaran2013mathcal][@ambikasaran2014fast]. The goal of ``HODLRlib`` is to serve as a reference implementation and is designed to have a simple interface with minimum number of dependencies. Our benchmarks demonstrate the large order of speedup that our library offers when compared to naive matrix factorizations.

Our code makes use of shared-memory parallelism through OpenMP. The solver is fairly generic and can handle matrices not necessarily arising out of kernel functions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the original articles.

``HODLRlib`` is designed such that the matrix corresponding to the linear system to be solved is abstracted through the ``HODLR_Matrix`` object, which needs to have the function ``getMatrixEntry``. This function takes in the arguments as the index in the matrix and returns the particular entry, which facilitates the reduction in storage costs, since only a few entries of the low-rank sub-blocks are needed to reconstruct these sub-blocks. This instance of ``HODLR_Matrix`` is then passed to the ``HODLR_Tree`` class, whose methods are used for the various matrix operations.

The current release has the following capabilities:

- MatVecs: Obtains $A x$ at a cost of $\mathcal{O}\left(N\log{N}\right)$
- Factorization: Factors the matrix $A$ into the desired form at a cost of $\mathcal{O}\left(N\log^2\left(N\right)\right)$
- Cholesky-like Symmetric Factorization: Obtains $A = W W^T$ at a cost of $\mathcal{O}\left(N\log^2\left(N\right)\right)$
- Solve: Solves linear systems $A x = b$ at an additional cost of $\mathcal{O}\left(N\log\left(N\right)\right)$
- Determinant Computation: Additional Cost of $\mathcal{O}\left(N\log{N} \right)$

``HODLRlib`` is released under the MPL2 license, and the source code is available at <https://github.com/sivaramambikasaran/HODLR>

# Acknowledgements

The authors were supported by the Department of Atomic Energy (DAE), India and Department of Science and Technology (DST), India.

# References
