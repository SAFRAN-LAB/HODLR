---
title: 'HODLRlib: An Library for HODLR Matrices'
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
date: 6th January 2019
bibliography: paper.bib
---

# Summary

``HODLRlib`` is a flexible library for working with matrices that have a Hierarchical Off-Diagonal Low-Rank (HODLR) structure. The current version performs matrix operations like matrix-vector products, solving and determinant computation in near-linear complexity, and has been demonstrated through our benchmarks. A key motivation for ``HODLRlib`` is to provide a high-performance, easy to use library for working with matrices that resemble a HODLR structure. 

``HODLRlib`` is an optimized implementation of the ideas illustrated in these articles[@ambikasaran2013mathcal][@ambikasaran2014fast]. The goal of ``HODLRlib`` is to serve as a reference implementation and is designed to have a simple interface with a minimum number of dependencies. Our benchmarks demonstrate the large order of magnitude speedup that our library offers when compared to naive matrix factorization.

Our code makes use of shared-memory parallelism through OpenMP. The solver has also been extended to matrices not necessarily arising out of kernels and also to higher dimensions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the original articles.

``HODLRlib`` is designed to be solved such that the matrix to be solved is abstracted through the ``HODLR_Matrix`` object which needs to have the function ``getMatrixEntry`` which takes in the arguments as the index in the matrix and returns the particular entry. This allows the reduction in storage costs. This instance of ``HODLR_Matrix`` is then passed to the ``HODLR_Tree`` class whose methods are used for the various matrix operations.

The current release has the following capabilities:

- Fast MatVecs: Obtains $A x$ at a cost of $\mathcal{O}\left(N\log{N}\right)$
- Fast Symmetric Factorization: Obtains $A = W W^T$ at a cost of $\mathcal{O}\left(N\log^2\left(N\right)\right)$
- Fast Solve: Solves linear systems $A x = b$ in $\mathcal{O}\left(N\log^2\left(N\right)\right)$
- Fast Determinant Computation: Additional Cost of $\mathcal{O}\left(N\log{N} \right)$

``HODLRlib`` is released under the MPL2 license, and the source code is available at <https://github.com/sivaramambikasaran/HODLR>

# Acknowledgements

# References
