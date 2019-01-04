********
Tutorial
********

For the sake of this tutorial, we are going to be using the ``tutorial.cpp`` file that is listed under ``examples/`` since it demonstrates all the features of this library. For the most part, comments in the file demonstrate intended functionality. However, we go over the main functions that may be of interest to a user in this page. 

**NOTE**: It is assumed that you have already completed the installation process of getting the dependencies.

Setting Parameters in CMakeLists.txt
------------------------------------

There are some variables that need to be set by the user at the top of the ``CMakeLists.txt`` file:

- ``INPUT_FILE``: This is the input ``.cpp`` file that needs to be compiled. For this tutorial, its going to be set to ``examples/tutorial.cpp``.
- ``OUTPUT_EXECUTABLE``: This is the name that the final build executable is given. Here we are just going to set is as ``tutorial``.
- ``DTYPE``: Datatype that is used in all the computations. Can be set to ``float``, ``double``, ``complex32`` and ``complex64``. We are just going to be using ``double`` for this tutorial.

Creating a Derived Class of ``HODLR_Matrix``:
---------------------------------------------

The matrix that needs to be solved for is abstracted through this derived class of ``HODLR_Matrix``. For the sake of the tutorial, we are calling this derived class ``Kernel``. The main method that needs to be set for this class is ``getMatrixEntry`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix. For instance for the Hilbert matrix, this would be set as::

    dtype getMatrixEntry(int i, int j) 
    {
        return (1 / (i + j + 1));
    }

Note that here ``dtype`` is set during compilation depending on ``DTYPE`` that was set in ``CMakeLists.txt``.

In this tutorial, we have initialized a random set of points in :math:`(-1, 1)` which are then sorted to obtain a coordinate vector :math:`x`. Using this we compute the distance between the :math:`i^{\mathrm{th}}` point and :math:`j^{\mathrm{th}}` in :math:`x` to obtain :math:`R` which can then be used in different kernel functions. 

Creating the Instance of ``HODLR_Tree``:
----------------------------------------

The main operations of this library are carried out through the ``HODLR_Tree`` class. The parameters that are taken for the constructor are the number of levels, tolerance for approximation and the earlier created instance of ``Kernel``::
    
    Kernel* K     = new Kernel(N, dim);
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, K);

We will now proceed to demonstrate the individual methods available under this class.

``assembleTree``
^^^^^^^^^^^^^^^^

We proceed to call the ``assembleTree`` method. This obtains the complete matrix for the leaf levels and the low-rank approximation for the off-diagonal blocks. Here we have used mentioned the fact that the matrix that we are constructing is both symmetric and positive-definite. Note than when we mention that the matrix is symmetric and positive-definite, the fast symmetric factorization method would be used. In all other cases the fast factorization method gets used::

    bool is_sym = true;
    bool is_pd = true;
    T->assembleTree(is_sym, is_pd);

``matmatProduct``
^^^^^^^^^^^^^^^^^

This function is used to obtain the matrix-matrix / matrix-vector product of the given matrix / vector :math:`x`, with the matrix that is abstracted by the instance of ``Kernel``::
    
    b = T->matmatProduct(x);

``factorize``
^^^^^^^^^^^^^

Depends upon whether we intend to perform fast factorization, or fast symmetric factorization:

- **Fast Factorization** - This function performs the factorizations such that the matrix is obtained as :math:`K = K_{\kappa} K_{\kappa-1} ... K_{1} K_{0}` where :math:`K_i` are block diagonal matrices with :math:`\kappa` being the number of levels considered. 

**Fast Symmetric Factorization** - This function performs the factorizations such that the matrix is obtained as :math:`K = K_{\kappa} K_{\kappa-1} ... K_{1} K_{0} K_{0}^T K_{1}^T ... K_{\kappa-1}^T K_{\kappa}^T` where :math:`K_i` are block diagonal matrices with :math:`\kappa` being the number of levels considered. 

For more details on this factorization refer to the articles `[1] <https://link.springer.com/article/10.1007/s10915-013-9714-z>`_ `[2] <https://arxiv.org/abs/1405.0223>`_

``solve``
^^^^^^^^^

Applies the inverse of the matrix(abstracted by the ``Kernel`` object) on the given vector :math:`x`. This must be called only after ``factorize`` has been called::

    x = T->solve(b);

``logDeterminant``
^^^^^^^^^^^^^^^^^^

Returns the log of the determinant of the matrix that has been described through the ``Kernel`` object::

    dtype log_det = T->logDeterminant();

``symmetricFactorProduct``
^^^^^^^^^^^^^^^^^^^^^^^^^^

If the matrix described through the ``Kernel`` object is a covariance matrix :math:`Q` it can be expressed as :math:`Q=W W^T`. If we create a random normal vector :math:`x` i.e :math:`\mathcal{N}(\mu = 0, \sigma = 1)`, then the random vector :math:`y` with covariance matrix :math:`Q` is given by :math:`y = W x`::

    y = T->symmetricFactorProduct(x);
