      program main

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Example program where fortran (fortran 90 here, using a feature
c     of fortran 2003) calls wrappers for HODLR C++ routines.
c
c     See: fort_hodlr_wrappers.cpp and fort_hodlr_wrappers.hpp for 
c     definitions of the c++ routines and a C wrapper for the C++ 
c     that prevents the linker from mangling the name of the 
c     subroutines.
c
c     To define the system matrix, edit the 
c     Matrix_Entry_Routine.cpp file.
c
c     Use accompanying makefile fort_hodlr_wrappers_dr.mk to
c     compile the code.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     The iso_c_binding is technically a fortran 2003 feature.
c     This works at least with gfortran (and many fortran 95
c     compilers)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use, intrinsic :: iso_c_binding
      
      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Interface tells fortran what to expect in terms of data types.
c     These "C" routines are defined in the fort_hodlr_wrappers.cpp
c     and fort_hodlr_wrappers.hpp files.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      interface
         subroutine initialize_matrix_wrap(A, kernel, N, nDim, nLeaf,
     1        xyz, xyztarg, work, diag, eps) bind(c)
         use iso_c_binding
         type(c_ptr) A
         type(c_ptr) kernel
         
         real( c_double ) xyz(*)
         real( c_double ) xyztarg(*)
         real( c_double ) work(*)
         real( c_double ) diag(*)
         
         integer( c_int ), value :: N
         integer( c_int ), value :: nDim
         integer( c_int ), value :: nLeaf

         real ( c_double ), value :: eps
         end subroutine
      end interface

      interface
         subroutine matrix_multiply_wrap(A, x, b, N, 
     1        nRow, nCol) bind(c)
         use iso_c_binding

         type(c_ptr) A
         
         real( c_double ) x(*)
         real( c_double ) b(*)
         
         integer( c_int ), value :: N
         integer( c_int ), value :: nRow
         integer( c_int ), value :: nCol

         end subroutine
      end interface

      interface
         subroutine matrix_factor_wrap(A) bind(c)
         use iso_c_binding

         type(c_ptr) A

         end subroutine
      end interface

      interface
         subroutine matrix_solve_wrap(A, x, b, N, 
     1        nRow, nCol) bind(c)
         use iso_c_binding

         type(c_ptr) A
         
         real( c_double ) x(*)
         real( c_double ) b(*)
         
         integer( c_int ), value :: N
         integer( c_int ), value :: nRow
         integer( c_int ), value :: nCol

         end subroutine
      end interface

      interface
         subroutine matrix_determinant_wrap(A, determinant) bind(c)
         use iso_c_binding

         type(c_ptr) A
         
         real( c_double ) determinant

         end subroutine
      end interface


      integer ( c_int ), parameter :: N = 100 000

      integer i

      real( c_double ) xyz(N)
      real( c_double ) xyztarg(N)
      real( c_double ) diag(N)
      real( c_double ) work(N)

      real ( c_double ) x(N)
      real ( c_double ) b(N)
      real ( c_double ) bResult(N)
      real ( c_double ) c(N)         
      real ( c_double ) cResult(N)         

      real ( c_double ) determinant

      type( c_ptr ) A
      type( c_ptr ) Atarg
      type( c_ptr ) kernel

      integer ( c_int ) nDim, nLeaf, nRow, nCol

      real( c_double ) eps

      real *8 t1, t2, res, resTotal, eInterp, sumInterp

c     set precision, number of dimensions, nLeaf, etc

      eps = 1.0d-15
      nDim = 1
      nLeaf = 100
      nRow = N
      nCol = 1
      determinant = 0.0d0

c     throw some values into these arrays
c     the xyz coordinates here are 1D and are already in sorted order

c     this performs the fit of y = sin(x) on [0,1] using the RBF
c     defined in Matrix_Entry_Routine.cpp at equispaced points
c     xyz and evaluates at the midpoints xyztarg

c     work is used to store a shape parameter for the RBF
      do i = 1,N
         diag(i) = 0.0d0
         xyz(i) = i*1.0d0/N
         xyztarg(i) = xyz(i)+0.5d0/N
         work(i) = N/100.0d0
         x(i) = i*0.5d0
         b(i) = sin(xyz(i))
         c(i) = sin(xyztarg(i))
      enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     call the c++ routines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write(*,*) 'Initializing matrix A ...'

      t1 = second()

      call initialize_matrix_wrap(A, kernel, N, nDim, nLeaf,
     1        xyz, xyz, work, diag, eps)

      t2 = second()

      write(*,*) 'time = ', t2-t1

      write(*,*) 'Initializing matrix Atarg ...'

      t1 = second()

      call initialize_matrix_wrap(Atarg, kernel, N, nDim, nLeaf,
     1        xyz, xyztarg, work, diag, eps)

      t2 = second()

      write(*,*) 'time = ', t2-t1


      write(*,*) 'Factoring Matrix A ...'

      t1 = second()

      call matrix_factor_wrap(A)

      t2 = second()

      write(*,*) 'time = ', t2-t1


      write(*,*) 'Factoring Matrix Atarg ...'

      t1 = second()

      call matrix_factor_wrap(Atarg)

      t2 = second()

      write(*,*) 'time = ', t2-t1


      write(*,*) 'Solving system x = A^(-1)b ...'

      t1 = second()

      call matrix_solve_wrap(A, x, b, N, nRow, nCol)

      t2 = second()

      write(*,*) 'time = ', t2-t1


      write(*,*) 'Calculating value of interpolant at targets ...'

      t1 = second()

      call matrix_multiply_wrap(Atarg, x, cResult, N, nRow, nCol)

      t2 = second()

      write(*,*) 'time = ', t2-t1


      write(*,*) 'and at interpolation points ...'

      t1 = second()

      call matrix_multiply_wrap(A, x, bResult, N, nRow, nCol)

      t2 = second()

      write(*,*) 'time = ', t2-t1


      write(*,*) 'Calculating log-determinant of A ...'

      t1 = second()

      call matrix_determinant_wrap(A, determinant)

      t2 = second()

      write(*,*) 'time = ', t2-t1

      write(*,*) 'determinant = ', determinant

      write(*,*) 'Calculating errors ... '

      res = 0.0d0
      resTotal = 0.0d0
      eInterp = 0.0d0
      sumInterp = 0.0d0

      do i = 1,N
         res = res + (b(i)-bResult(i))**2
         resTotal = resTotal + b(i)**2
         eInterp = eInterp + (c(i)-cResult(i))**2
         sumInterp = sumInterp + c(i)**2
      enddo

      write(*,*) 'Residual error ', dsqrt(res)

      write(*,*) 'Relative residual error ', dsqrt(res/resTotal)
      
      write(*,*) 'Interpolation error at targets ', dsqrt(eInterp)

      write(*,*) 'Relative interpolation error at targets ', 
     1     dsqrt(eInterp/sumInterp)

      stop
      end
