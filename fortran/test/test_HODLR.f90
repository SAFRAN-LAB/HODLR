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
    ! Log-Det exact for this matrix with N = 1000:
    real(c_double) :: log_det_exact = 4598.3043670568

    ! Method used for low-rank matrix decomposition:
    character(len = 20) :: factorization_method = "rookPivoting"

    ! Creating the kernel object:
    call initialize_kernel_object(N, kernel)

    ! Getting the matrix encoded by the kernel object:
    call get_matrix(kernel, 0, 0, N, N, flattened_matrix)
    matrix = reshape(flattened_matrix, (/ N, N /))

    ! Building the HODLR object:
    call initialize_hodlr(N, M, eps, hodlr)
    call assemble(hodlr, kernel, "rookPivoting", is_sym, is_pd)

    ! Populating the vector x with random elements:
    do i = 1, N
        call random_number(x(i)) 
    end do
    
    ! Performing matrix multiplication using HODLR:
    call matmat_product(hodlr, x, b)
    ! Performing the matrix multiplication directly:
    b_exact = matmul(matrix, x)
    ! Checking the error of matmul:
    if((sum(abs(b_exact - b)) / sum(abs(b_exact)) > N * eps)) then
        print *, "POSSIBLE BUG!!: Large error in MatMul"
        call exit(1)
    end if

    ! Factorizing the elements of the tree:
    call factorize(hodlr)

    ! Solve for x in the equation Ax = b:
    call solve(hodlr, b_exact, x_approx)

    ! Checking the error of solve:
    if((sum((x_approx - x)**2) / sum(x**2) > N * eps)) then
        print *, "POSSIBLE BUG!!: Large error in Solve(Forward)"
        call exit(1)
    end if

    if((sum((matmul(matrix, x) - b_exact)**2) / sum(b_exact**2) > N * eps)) then
        print *, "POSSIBLE BUG!!: Large error in Solve(Backward)"
        call exit(1)
    end if

    ! Getting the determinant of the matrix encoded:
    call logdeterminant(hodlr, log_det)
    if(abs(log_det - log_det_exact) / log_det_exact > 1e-7) then
        print *, "POSSIBLE BUG!!: Large error in log determinant"
        call exit(1)
    end if

    ! We set y = W^T x
    call symmetric_factor_transpose_product(hodlr, x, y)
    ! b = W y = W W^T x = B * x
    call symmetric_factor_product(hodlr, y, b)
    if(sum((b_exact - b)**2) / sum(b_exact**2) > N * eps) then
        print *, "POSSIBLE BUG!!: Large error in symmetric factor products"
        call exit(1)
    end if

    ! Getting the matrix encoded by the kernel object:
    call get_symmetric_factor(hodlr, flattened_W_matrix)
    W_matrix = reshape(flattened_W_matrix, (/ N, N /))
    if(sum(abs(matrix - matmul(W_matrix, transpose(W_matrix)))) / sum(matrix) > N * eps) then
        print *, "POSSIBLE BUG!!: Large error in obtaining symmetric factor"
        call exit(1)
    end if

    print *, "Reached End of Test File Successfully! All functions work as intended!"
end
