module hodlr_mod
    
    interface
        subroutine initialize_kernel_object(kernel, N, dim) &
        bind(c, name = "initialize_kernel_object_c")
        
            use iso_c_binding

            type(c_ptr) :: kernel

            integer(c_int), value :: N
            integer(c_int), value :: dim

        end subroutine
    end interface

    interface
        subroutine get_matrix(matrix, kernel, row_start, col_start, row_end, col_end) &
        bind(c, name = "get_matrix_c")
        
            use iso_c_binding

            real(c_double) :: matrix(*)
            type(c_ptr) :: kernel

            integer(c_int), value :: row_start
            integer(c_int), value :: col_start
            integer(c_int), value :: row_end
            integer(c_int), value :: col_end
       end subroutine
    end interface

    interface
        subroutine initialize_matrix_factorizer(factorizer, kernel, factorization_method) &
        bind(c, name = "initialize_matrix_factorizer_c")
        
            use iso_c_binding

            type(c_ptr) :: factorizer
            type(c_ptr) :: kernel

            character(c_char) :: factorization_method(*)
        end subroutine
    end interface

    interface
        subroutine get_factorization(factorizer, l, r, eps) &
        bind(c, name = "get_factorization_c")

            use iso_c_binding

            type(c_ptr) :: factorizer

            real(c_double) :: l(*)
            real(c_double) :: r(*)
            real(c_double), value :: eps

        end subroutine
    end interface

    interface
        subroutine initialize_hodlr_tree(tree, n_levels, eps, factorizer) &
        bind(c, name = "initialize_hodlr_tree_c")
            
            use iso_c_binding

            type(c_ptr) :: tree

            integer(c_int), value :: n_levels
            real(c_double), value :: eps

            type(c_ptr) :: factorizer

        end subroutine
    end interface

    interface
        subroutine assemble_tree(tree, is_sym, is_pd) &
        bind(c, name = "assemble_tree_c")

            use iso_c_binding

            type(c_ptr) ::tree

            logical(c_bool), value :: is_sym
            logical(c_bool), value :: is_pd

        end subroutine
    end interface

    interface
        subroutine matmat_product(tree, x, b) &
        bind(c, name = "matmat_product_c")

            use iso_c_binding

            type(c_ptr) ::tree
            real(c_double) :: x(*)
            real(c_double) :: b(*)

        end subroutine
    end interface

    interface
        subroutine factorize(tree) &
        bind(c, name = "factorize_c")

            use iso_c_binding
            type(c_ptr) ::tree

        end subroutine
    end interface

    interface
        subroutine solve(tree, b, x) &
        bind(c, name = "solve_c")

            use iso_c_binding

            type(c_ptr) :: tree
            real(c_double) :: b(*)
            real(c_double) :: x(*)

        end subroutine
    end interface

    interface
        subroutine logdeterminant(tree, log_det) &
        bind(c, name = "logdeterminant_c")

            use iso_c_binding

            type(c_ptr) :: tree
            real(c_double) :: log_det

        end subroutine
    end interface

end module hodlr_mod
