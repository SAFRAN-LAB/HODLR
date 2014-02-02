Python bindings to HODLR solver
===============================

This is a light Python wrapper around the HODLR solver.

Installation
------------

In the ``python`` directory (the location of this file), run

::

    python setup.py install

or, if you've installed Eigen in a custom location, run

::

    python setup.py build_ext --include-dirs=/path/to/eigen3
    python setup.py install

Usage
-----

This module exposes one class ``HODLR`` to perform all the computations. The
basic usage will look like:

::

    import hodlr
    import numpy as np

    def matrix(i, j):
        return np.exp(-0.5 * (i - j) ** 2)

    # Set up the system.
    N = 500
    diag = np.random.randn(N)
    solver = hodlr.HODLR(matrix, diag)

    # Build the true solution for testing purposes.
    x = np.random.randn(N)
    b = solver.matrix_product(x)

    # Solve the system.
    xsol = solver.solve(b)
    print("L1 error: {0}".format(np.mean(np.abs(x - xsol))))

    # Get the log-determinant.
    print("Log-determinant: {0}".format(solver.logdet()))

Take a look at ``demo.py`` for a more detailed example.

Authors
-------

The solver is developed by Sivaram Ambikasaran and these Python bindings are
developed by Daniel Foreman-Mackey.
