module libkernel
    use iso_c_binding

    private
    public :: kernel

    include "kernel_cdef.f90"

    type kernel
        private
        type(c_ptr) :: ptr 
    contains

        final :: delete_kernel ! Destructor
        procedure :: get_matrix => foo_baz

    end type

    ! This function will act as the constructor for kernel type
    interface kernel
        procedure create_kernel
    end interface

! Implementation of the functions. We just wrap the C function here.
contains 
    function create_kernel(N)
        implicit none
        type(kernel) :: create_kernel
        integer, intent(in) :: N
        create_kernel%ptr = create_kernel(a, b)
    end function

    subroutine delete_foo(this)
        implicit none
        type(foo) :: this
        call delete_foo_c(this%ptr)
    end subroutine

    double precision function foo_baz(this, c)
        implicit none
        class(foo), intent(in) :: this
        double precision, intent(in) :: c
        foo_baz = foo_baz_c(this%ptr, c)
    end function

    subroutine foo_speaker(str)
        implicit none
        character(len=*), intent(in) :: str
        character(len=1, kind=C_CHAR) :: c_str(len_trim(str) + 1)
        integer :: N, i

        ! Converting Fortran string to C string
        N = len_trim(str)
        do i = 1, N
            c_str(i) = str(i:i)
        end do
        c_str(N + 1) = C_NULL_CHAR

        call foo_speaker_c(c_str)
    end subroutine
end module
