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

    ! Size of the matrix:
    integer, parameter :: N = 100
    ! Rank of the factorization performed:
    integer :: rank
    ! Number of digits of tolerance required:
    real(c_double), parameter :: eps = 1.0d-14
    ! Used to store the complete matrix:
    real(c_double) :: flattened_matrix(N * N)
    real(c_double) :: matrix(N, N)

    ! The following arrays are used in getting the direct factorization of the array:
    real(c_double) :: l_flat(N * N)
    real(c_double) :: r_flat(N * N)
    real(c_double), dimension (:,:), allocatable :: l
    real(c_double), dimension (:,:), allocatable :: r
    
    integer :: i, j

    ! Method used for low-rank matrix decomposition:
    character(len = 20) :: factorization_method = "SVD"

    ! Creating the kernel object:
    call initialize_kernel_object(N, kernel)

    ! Getting the matrix encoded by the kernel object:
    call get_matrix(kernel, 0, 0, N, N, flattened_matrix)
    matrix = reshape(flattened_matrix, (/ N, N /))

    ! Building the matrix factorizer object:
    call initialize_matrix_factorizer(kernel, factorization_method, factorizer)

    ! ! Example of directly getting the factorization:
    call get_factorization(factorizer, eps, l_flat, r_flat, rank)

    allocate(l(N, rank))
    allocate(r(N, rank))

    ! Reshaping the flattened matrices:
    l = reshape(l_flat(1:N * rank), (/ N, rank /))
    r = reshape(r_flat(1:N * rank), (/ N, rank /))

    ! Printing the error:
    print *, "Error in Factorization:", SUM(ABS(matrix - matmul(l, transpose(r))))

    deallocate(l)
    deallocate(r)

end
