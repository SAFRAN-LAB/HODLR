program main
   use iso_c_binding
   implicit none

   interface
      subroutine initialize_kernel_object(kernel, N, dim) bind(c)
         use iso_c_binding

         type(c_ptr), intent(inout) :: kernel

         integer(c_int), value :: N
         integer(c_int), value :: dim

      end subroutine
   end interface

   interface
      subroutine get_matrix(matrix, kernel, row_start, col_start, row_end, col_end) bind(c)
         use iso_c_binding

         real(c_double) :: matrix(*)
         type(c_ptr) :: kernel

         integer(c_int), value :: row_start
         integer(c_int), value :: col_start
         integer(c_int), value :: row_end
         integer(c_int), value :: col_end
      end subroutine
   end interface

   type(c_ptr) :: kernel
   type(c_ptr) :: factorizer
   type(c_ptr) :: tree

   integer, parameter :: N = 5
   integer, parameter :: dim = 1

   ! Used to store the complete matrix:
   real(c_double) :: flattened_matrix(N * N)
   real(c_double) :: matrix(N, N)

   integer :: i, j

   ! Now, we call the C++ routines:
   print *, "Initializing the Kernel Object"
   call initialize_kernel_object(kernel, N, dim)
   print *, "Successfully created Kernel Object"
   print *, "Building the matrix..."
   call get_matrix(flattened_matrix, kernel, 0, 0, N, N)
   print *, "Successfully built the matrix"
   matrix = reshape(flattened_matrix, (/ N, N /))
   print *, "Printing the Matrix:"
   do i = 1, 5
       print *, (matrix(i, j), j = 1, 5)
   end do

end
