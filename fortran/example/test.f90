program test
    
    use libkernel
    implicit none
    type(kernel) :: k

    ! Size of the matrix considered
    integer, parameter :: N = 5
    ! Used to recover the whole matrix
    double precision :: a(N, N)
    ! Counters used to loop over elements of the data matrix:
    integer :: i, j

    double precision :: vector_x(N)

    ! Random number:
    double precision :: rand

    ! Create an object of type foo
    k = kernel(N)
    print *, k%getmatrixentry(0, 0)


    call k%getvectorx(p)

    ! call k%getmatrix(0, 0, N, N, a)

! The destructor should be called automatically here, but this is not yet
! implemented in gfortran. So let's do it manually.
#ifdef __GNUC__
    call k%delete
#endif

end program
