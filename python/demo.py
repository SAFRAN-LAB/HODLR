#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import time
import hodlr
import numpy as np

def matrix(i, j):
    return np.exp(-0.5 * (i - j) ** 2)

# Set up the system.
N = 5000
diag = np.random.randn(N)
print("Setting up...")
strt = time.time()
solver = hodlr.HODLR(matrix, diag)
print("Took {0} seconds".format(time.time()-strt))

# Build the true solution for testing purposes.
x = np.random.randn(N)
print("Computing matrix-matrix product...")
strt = time.time()
b = solver.matrix_product(x)
print("Took {0} seconds".format(time.time()-strt))

# Solve the system.
print("Solving...")
strt = time.time()
xsol = solver.solve(b)
print("Took {0} seconds".format(time.time()-strt))
print("L1 error: {0}".format(np.mean(np.abs(x - xsol))))

# Get the log-determinant.
print("Finding log-determinant...")
strt = time.time()
logdet = solver.logdet()
print("Took {0} seconds".format(time.time()-strt))
print("Value: {0}".format(logdet))
