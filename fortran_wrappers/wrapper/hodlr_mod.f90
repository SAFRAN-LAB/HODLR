module hodlr_mod
    
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

            real(c_double), intent(inout) :: matrix(*)
            type(c_ptr), intent(in) :: kernel

            integer(c_int), value :: row_start
            integer(c_int), value :: col_start
            integer(c_int), value :: row_end
            integer(c_int), value :: col_end
       end subroutine
    end interface

    interface
        subroutine initialize_matrix_factorizer(factorizer, kernel, factorization_method) bind(c)
            use iso_c_binding

            type(c_ptr), intent(inout) :: factorizer
            type(c_ptr), intent(in) :: kernel

            character(c_char), intent(in) :: factorization_method(*)
        end subroutine
    end interface

    interface
        subroutine get_factorization(factorizer, l, r, eps) bind(c)
            use iso_c_binding

            type(c_ptr), intent(in) :: factorizer

            real(c_double), intent(inout) :: l(*)
            real(c_double), intent(inout) :: r(*)
            real(c_double), intent(in) :: eps

        end subroutine
    end interface

    ! interface
    !     subroutine initialize_hodlr_tree(tree, n_levels, eps, factorizer) bind(c)
    !         use iso_c_binding

    !         type(c_ptr), intent(inout) :: tree

    !         integer(c_int), intent(in) :: n_levels
    !         real(c_double), intent(in) :: eps

    !         type(c_ptr), intent(in) :: factorizer

    !     end subroutine
    ! end interface

end module hodlr_mod
