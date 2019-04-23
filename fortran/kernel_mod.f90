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
        ! procedure :: getmatrix => get_matrix

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

    function get_vector_x(this)
        implicit none
        class(kernel), intent(in) :: this
        double precision, pointer :: get_vector_x
        
        get_vector_x = get_vector_x_c(this%ptr)
    end function

    ! subroutine get_matrix(this, row_start, col_start, row_end, col_end, a)
    !     implicit none
    !     class(kernel), intent(in) :: this
    !     integer, intent(in) :: row_start
    !     integer, intent(in) :: col_start
    !     integer, intent(in) :: row_end
    !     integer, intent(in) :: col_end

    !     double precision :: a(row_end - row_start, col_end - col_start)
    !     call get_matrix_c(this%ptr, row_start, col_start, row_end, col_end, a)
    ! end subroutine

end module
