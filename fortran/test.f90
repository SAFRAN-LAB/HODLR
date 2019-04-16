program test
    
    use libkernel
    implicit none
    type(kernel) :: k

    ! Create an object of type foo
    k = kernel(10)

    print *, k%getmatrixentry(0, 0)

! The destructor should be called automatically here, but this is not yet
! implemented in gfortran. So let's do it manually.
#ifdef __GNUC__
    call k%delete
#endif

end program
