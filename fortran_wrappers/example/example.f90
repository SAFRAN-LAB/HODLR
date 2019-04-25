program main

   use iso_c_binding
   use hodlr_mod

   implicit none

   ! The following objects are used to refer to the C++ objects
   type(c_ptr) :: kernel
   type(c_ptr) :: factorizer
   type(c_ptr) :: tree

   ! Size of the matrix:
   integer, parameter :: N = 10000
   ! Size of the matrix at the leaf level:
   integer, parameter :: M = 200
   ! Dimensionality of the problem considered:
   integer, parameter :: dim = 1
   ! Number of levels in the tree:
   ! integer :: n_levels = LOG(REAL(N) / REAL(M)) / LOG(2.)
   ! Tolerance for the approximation used:
   real(c_double), parameter :: eps = 1e-12

   ! Used to store the complete matrix:
   real(c_double) :: flattened_matrix(N * N)
   real(c_double) :: matrix(N, N)
   ! The following arrays are used in getting the direct factorization of the array:
   real(c_double) :: l_flat(N * N)
   real(c_double) :: r_flat(N * N)
   real(c_double) :: l(N, N)
   real(c_double) :: r(N, N)

   integer :: i, j

   ! Method used for low-rank matrix decomposition:
   character(len = 20) :: factorization_method = "rookPivoting"

   ! Creating the kernel object:
   print *, 'HERE'
   call initialize_kernel_object(kernel, N, dim)
   print *, 'HERE2'

   ! Getting the matrix encoded by the kernel object:
   call get_matrix(flattened_matrix, kernel, 0, 0, N, N)
   print *, 'HERE3'
   matrix = reshape(flattened_matrix, (/ N, N /))
   print *, 'HERE4'

   ! Printing the Matrix to observe the structure:
   ! call print_matrix(matrix, N)

   ! Building the matrix factorizer object:
   call initialize_matrix_factorizer(factorizer, kernel, factorization_method)
   print *, 'HERE5'

   ! Example of directly getting the factorization:
   ! call get_factorization(factorizer, l_flat, r_flat, eps)

   ! Reshaping the flattened matrices:
   ! l = reshape(l_flat, (/ N, N /))
   ! r = reshape(r_flat, (/ N, N /))

   ! Printing the Matrices and the error matrix:
   ! call print_matrix(l, N)
   ! call print_matrix(r, N)
   ! call print_matrix(matrix - matmul(l, r), N)

   ! Building the HODLR tree object:
   ! call initialize_hodlr_tree(tree, n_levels, eps, factorizer)

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
