! Function that returns the matrix values:
function get_matrix_entry(i, j) bind(c, name="get_matrix_entry")
    implicit none
    integer, intent(in) :: i
    integer, intent(in) :: j

    ! It needs to be ensured that this is the same as N set below:
    integer, parameter :: N = 1000
    double precision :: get_matrix_entry

    if(i == j) then
        get_matrix_entry = 100
    else
        get_matrix_entry = exp(-(DBLE(i)/N - DBLE(j)/N)**2)
    end if

end function

program main

    use iso_c_binding
    use hodlr_mod

    implicit none

    ! The following objects are used to refer to the C++ objects
    type(c_ptr) :: kernel
    type(c_ptr) :: hodlr

    ! Size of the matrix:
    integer, parameter :: N = 1000
    ! Size of the matrix at leaf level:
    integer, parameter :: M = 200
    ! Number of digits of tolerance required:
    real(c_double), parameter :: eps = 1.0d-15

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
    character(len = 20) :: factorization_method = "SVD"

    ! Creating the kernel object:
    call initialize_kernel_object(N, kernel)

    ! Getting the matrix encoded by the kernel object:
    call get_matrix(kernel, 0, 0, N, N, flattened_matrix)
    matrix = reshape(flattened_matrix, (/ N, N /))

    ! Printing the Matrix to observe the structure:
    ! print *, "Printing the Encoded Matrix A:"
    ! call print_matrix(matrix, N)

    ! Building the HODLR object:
    call initialize_hodlr(N, M, eps, kernel, "rookPivoting", is_sym, is_pd, hodlr)

    ! Populating the vector x with random elements:
    do i = 1, N
        call random_number(x(i)) 
    end do
    
    ! Performing matrix multiplication using HODLR:
    call matmat_product(hodlr, x, b)
    ! Performing the matrix multiplication directly:
    b_exact = matmul(matrix, x)
    ! Printing the error of matmul:
    print *, "Error in Matrix Multiplication:", (sum(abs(b_exact - b)) / sum(abs(b_exact)))

    ! Factorizing the elements of the tree:
    call factorize(hodlr)

    ! Solve for x in the equation Ax = b:
    call solve(hodlr, b_exact, x_approx)

    ! Printing the error of solve:
    print *, "Error in Solve(Forward):", sum((x_approx - x)**2) / sum(x**2)
    print *, "Error in Solve(Backward):", sum((matmul(matrix, x) - b_exact)**2) / sum(b_exact**2)

    ! Getting the determinant of the matrix encoded:
    call logdeterminant(hodlr, log_det)
    print *, "Log Determinant of Matrix:", log_det

    ! We set y = W^T x
    call symmetric_factor_transpose_product(hodlr, x, y)
    ! b = W y = W W^T x = B * x
    call symmetric_factor_product(hodlr, y, b)
    print *, "Error in Matrix Multiplication(using Symmetric Factor Products):", sum((b_exact - b)**2) / sum(b_exact**2)

    ! Getting the matrix encoded by the kernel object:
    call get_symmetric_factor(hodlr, flattened_W_matrix)
    W_matrix = reshape(flattened_W_matrix, (/ N, N /))
    print *, "Error in recovering A from W", sum(abs(matrix - matmul(W_matrix, transpose(W_matrix)))) / sum(matrix)
    
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
