! Function that returns the matrix values:
function get_matrix_entry(i, j) bind(c, name="get_matrix_entry")
    implicit none
    integer, intent(in) :: i
    integer, intent(in) :: j

    double precision :: get_matrix_entry
    get_matrix_entry = 1 / DFLOAT(i + j + 1)
end function

program main

    use iso_c_binding
    use hodlr_mod

    implicit none

    ! The following objects are used to refer to the C++ objects
    type(c_ptr) :: kernel
    type(c_ptr) :: factorizer
    type(c_ptr) :: tree

    ! Size of the matrix:
    integer, parameter :: N = 5
    ! Size of the matrix as 5 at leaf level:
    integer, parameter :: M = 5
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

    ! Variables used to store the symmetric factor matrix:
    real(c_double) :: flattened_W_matrix(N * N)
    real(c_double) :: W_matrix(N, N)
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
    real(c_double) :: y(N)

    ! Used to store the log determinant of the matrix:
    real(c_double) :: log_det

    ! Method used for low-rank matrix decomposition:
    character(len = 20) :: factorization_method = "rookPivoting"

    ! Creating the kernel object:
    call initialize_kernel_object(N, kernel)

    ! Getting the matrix encoded by the kernel object:
    call get_matrix(kernel, 0, 0, N, N, flattened_matrix)
    matrix = reshape(flattened_matrix, (/ N, N /))

    ! Printing the Matrix to observe the structure:
    print *, "Printing the Encoded Matrix A:"
    call print_matrix(matrix, N)

    ! Building the matrix factorizer object:
    call initialize_matrix_factorizer(kernel, factorization_method, factorizer)

    ! ! Example of directly getting the factorization:
    call get_factorization(factorizer, eps, l_flat, r_flat)

    ! Reshaping the flattened matrices:
    l = reshape(l_flat, (/ N, N /))
    r = reshape(r_flat, (/ N, N /))

    ! Printing the Matrices and the error matrix:
    print *, "Printing Matrix L:"
    call print_matrix(l, N)
    print *, "Printing Matrix R:"
    call print_matrix(r, N)
    print *, "Printing Error Matrix (A - L * R):"
    call print_matrix(matrix - matmul(l, r), N)

    ! Building the HODLR tree object:
    call initialize_hodlr_tree(n_levels, eps, factorizer, tree)
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
    print *, "Log Determinant of Matrix:", log_det

    ! We set y = W^T x
    call symmetric_factor_transpose_product(tree, x, y)
    ! b = W y = W W^T x = B * x
    call symmetric_factor_product(tree, y, b)
    print *, "Error in Matrix Multiplication(using Symmetric Factor Products):", sum(abs(b_exact - b))

    ! Getting the matrix encoded by the kernel object:
    call get_symmetric_factor(tree, flattened_W_matrix)
    W_matrix = reshape(flattened_W_matrix, (/ N, N /))
    print *, "Error in recovering A from W", sum(matrix - matmul(W_matrix, transpose(W_matrix)))
    
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
