program main

    use iso_c_binding
    use hodlr_mod

    implicit none

    ! The following objects are used to refer to the C++ objects
    type(c_ptr) :: kernel
    type(c_ptr) :: factorizer
    type(c_ptr) :: tree

    ! Size of the matrix:
    integer, parameter :: N = 1000
    ! Size of the matrix at the leaf level:
    integer, parameter :: M = 200
    ! Dimensionality of the problem considered:
    integer, parameter :: dim = 1
    ! Number of digits of tolerance required:
    real(c_double), parameter :: eps = 1.0d-15
    ! Number of levels in the tree:
    integer :: n_levels = LOG(REAL(N) / REAL(M)) / LOG(2.)

    ! Whether matrix encoded is symmetric:
    logical(c_bool) :: is_sym = .true.
    ! Whether matrix encoded is positive-definite:
    logical(c_bool) :: is_pd = .true.

    ! Used to store the complete matrix:
    real(c_double) :: flattened_matrix(N * N)
    real(c_double) :: matrix(N, N)
    ! The following arrays are used in getting the direct factorization of the array:
    real(c_double) :: l_flat(N * N)
    real(c_double) :: r_flat(N * N)
    real(c_double) :: l(N, N)
    real(c_double) :: r(N, N)

    integer :: i, j

    ! Creating random vectors used on the RHS and LHS for multiplication and solving:
    real(c_double) :: x(N)
    real(c_double) :: b(N)
    real(c_double) :: b_exact(N)
    real(c_double) :: x_approx(N)

    ! Used to store the log determinant of the matrix:
    real(c_double) :: log_det

    ! Method used for low-rank matrix decomposition:
    character(len = 20) :: factorization_method = "rookPivoting"

    ! Creating the kernel object:
    call initialize_kernel_object(kernel, N, dim)

    ! Getting the matrix encoded by the kernel object:
    call get_matrix(flattened_matrix, kernel, 0, 0, N, N)
    matrix = reshape(flattened_matrix, (/ N, N /))

    ! Printing the Matrix to observe the structure:
    ! call print_matrix(matrix, N)

    ! Building the matrix factorizer object:
    call initialize_matrix_factorizer(factorizer, kernel, factorization_method)

    ! Example of directly getting the factorization:
    ! call get_factorization(factorizer, l_flat, r_flat, eps)

    ! ! Reshaping the flattened matrices:
    ! l = reshape(l_flat, (/ N, N /))
    ! r = reshape(r_flat, (/ N, N /))

    ! ! Printing the Matrices and the error matrix:
    ! call print_matrix(l, N)
    ! call print_matrix(r, N)
    ! call print_matrix(matrix - matmul(l, r), N)

    ! Building the HODLR tree object:
    call initialize_hodlr_tree(tree, n_levels, eps, factorizer)
    ! Assembling the tree under awareness of positive definite and symmetric nature:
    call assemble_tree(tree, is_sym, is_pd)

    ! Populating the vector x with random elements:
    do i = 1, N
        call random_number(x(i)) 
    end do
    
    ! Performing matrix multiplication using HODLR:
    call matmat_product(tree, x, b)
    ! Performing the matrix multiplication directly:
    b_exact = matmul(matrix, x)
    ! Printing the error of matmul:
    print *, "Error in Matrix Multiplication:", sum(abs(b_exact - b))

    ! Factorizing the elements of the tree:
    call factorize(tree)

    ! Solve for x in the equation Ax = b:
    call solve(tree, b_exact, x_approx)

    ! Printing the error of solve:
    print *, "Error in Solve:", sum(abs(x_approx - x))

    ! Getting the determinant of the matrix encoded:
    call logdeterminant(tree, log_det)
    print *, "Determinant of Matrix:", log_det

    contains
        ! Use this subroutine to print and check the matrix:
        subroutine print_matrix(matrix, N)
            implicit none
            integer, intent(in) :: N
            real(c_double), intent(in) :: matrix(N, N)
            
            do i = 1, N
                print *, (matrix(i, j), j = 1, N)
            end do
        end subroutine
end
