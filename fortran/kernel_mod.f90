module libkernel

    use iso_c_binding

    private
    public :: kernel

    ! Yes, include is a keyword in Fortran !
    include "kernel_cdef.f90"

    ! We'll use a Fortan type to represent a C++ class here, in an opaque maner
    type kernel
        private
        type(c_ptr) :: ptr ! pointer to the Kernel class
    contains

! We can bind some functions to this type, allowing for a cleaner syntax.
#ifdef __GNUC__
        procedure :: delete => delete_kernel_polymorph ! Destructor for gfortran
#else
        final :: delete_kernel ! Destructor
#endif

        ! Function member
        procedure :: getmatrixentry => get_matrix_entry

    end type

    ! This function will act as the constructor for Kernel type
    interface kernel
        procedure create_kernel
    end interface

! Implementation of the functions. We just wrap the C function here
contains

    function create_kernel(N)
        implicit none
        type(kernel) :: create_kernel
        integer, intent(in) :: N
        create_kernel%ptr = create_kernel_c(N)
    end function

    subroutine delete_kernel(this)
        implicit none
        type(kernel) :: this
        call delete_kernel_c(this%ptr)
    end subroutine

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_kernel_polymorph(this)
        implicit none
        class(kernel) :: this
        call delete_kernel_c(this%ptr)
    end subroutine

    double precision function get_matrix_entry(this, i, j)
        implicit none
        class(kernel), intent(in) :: this
        integer, intent(in) :: i
        integer, intent(in) :: j

        get_matrix_entry = get_matrix_entry_c(this%ptr, i, j)
    end function

end module
