---
title: 'HODLRlib: A Library for Hierarchical Matrices'
tags:
  - Hierarchical Matrix
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

``HODLRlib`` is a flexible library for working with matrices that have a Hierarchical Off-Diagonal Low-Rank (HODLR) [@ambikasaran2013mathcal] structure. The current version performs matrix operations like matrix-vector products, solving linear systems, Cholesky-like symmetric factorization and determinant computation in almost linear complexity. 

Significant speedups are obtained with ``HODLRlib``, when dealing with matrices that have an underlying data sparsity (rank-structured). Most rank-structured matrices fall in the class of Hierarchical matrices: some examples are $\mathcal{H}$ [@hmatrix1999; @hmatrix2000], $\mathcal{H}^2$ [@h2matrix2002], HSS [@chandrasekaran2006fast], HODLR. These formats differ from each other in their matrix partitioning, low-rank sub-blocks and nested basis behaviour. A thorough comparison of these different formats has been provided in @ambikasaran2013fast.

Among these hierarchical matrices, HODLR matrices are applicable for a wide range of problems of practical interest [@kressner2017fast; @Hartmann2018; @ambikasaran2016fast]. In addition, they possess the simplest data-sparse representation, which inturn allows for easier implementation and obtaining great speed-ups.

``HODLRlib`` strives to provide a functionality for ``HODLR`` matrices similar to how libraries such as STRUMPACK [@ghysels2016efficient] and H2lib [@h2matrix2002] have provided for the HSS and $\mathcal{H}^2$ matrix formats respectively. Some unique features of our library are Cholesky-like symmetric factorization and determinant computation. While it is difficult to determine which matrix format works best, studies [@rouet_2015] in the past have revealed that it isn't a case of "one size fits all". Rather, it reveals that the optimum hierarchical matrix format is dependent on the problem and the size considered.

``HODLRlib`` is an optimized implementation of the ideas illustrated in @ambikasaran2013mathcal, @ambikasaran2016fast and @ambikasaran2014fast. The goal of ``HODLRlib`` is to serve as a high-performance, easy to use reference implementation for working with any HODLR matrix. We believe that apart from the features provided by the library, the simplicity and the extendability of the implementation is our USP. Our interfaces are kept simple with a minimum number of dependencies, with well-documented code to ensure that a potential user/developer can hit the ground running as soon as possible.

Our code makes use of shared-memory parallelism through OpenMP. The solver is fairly generic and can handle matrices not necessarily arising out of kernel functions. Further, the solver has been optimized and the running time of the solver is now massively (a few orders of magnitude) faster than the running times reported in the original articles.

``HODLRlib`` is designed such that the matrix corresponding to the linear system to be solved is abstracted through the ``HODLR_Matrix`` object, which needs to have the function ``getMatrixEntry``. This function takes in the arguments as the index in the matrix and returns the particular entry, which facilitates the reduction in storage costs, since only a few entries of the low-rank sub-blocks are needed to reconstruct these sub-blocks. This instance of ``HODLR_Matrix`` is then passed to the ``HODLR_Tree`` class, whose methods are used for the various matrix operations.

The current release has the following capabilities:

- MatVecs: Obtains $A x$ at a cost of $\mathcal{O}\left(N\log{N}\right)$
- Factorization: Factors the matrix $A$ into the desired form at a cost of $\mathcal{O}\left(N\log^2\left(N\right)\right)$
- Cholesky-like Symmetric Factorization: Obtains $A = W W^T$ at a cost of $\mathcal{O}\left(N\log^2\left(N\right)\right)$
- Solve: Solves linear systems $A x = b$ at an additional cost of $\mathcal{O}\left(N\log\left(N\right)\right)$
- Determinant Computation: Additional Cost of $\mathcal{O}\left(N\log{N} \right)$

``HODLRlib`` is released under the MPL2 license.

# Acknowledgements

The authors were supported by the Department of Atomic Energy (DAE), India and Department of Science and Technology (DST), India. The authors would also like to thanks Vaishnavi Gujjula and Michael Hartmann for their valuable comments, suggestions and code edits.

# References
