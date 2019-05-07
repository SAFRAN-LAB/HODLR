module hodlr_mod
    
    interface
        subroutine initialize_kernel_object(N, kernel) &
        bind(c, name = "initialize_kernel_object_c")
        
            use iso_c_binding

            integer(c_int), value :: N
            type(c_ptr) :: kernel

        end subroutine
    end interface

    interface
        subroutine get_matrix(kernel, row_start, col_start, row_end, col_end, matrix) &
        bind(c, name = "get_matrix_c")
        
            use iso_c_binding

            type(c_ptr) :: kernel

            integer(c_int), value :: row_start
            integer(c_int), value :: col_start
            integer(c_int), value :: row_end
            integer(c_int), value :: col_end
            
            real(c_double) :: matrix(*)
       
       end subroutine
    end interface

    interface
        subroutine initialize_lowrank(kernel, lowrank_method, lowrank) &
        bind(c, name = "initialize_lowrank_c")
        
            use iso_c_binding

            type(c_ptr) :: kernel
            character(c_char) :: lowrank_method(*)
            type(c_ptr) :: lowrank

        end subroutine
    end interface

    interface
        subroutine get_rank(lowrank, eps, rank) &
        bind(c, name = "get_rank_c")

            use iso_c_binding

            type(c_ptr) :: lowrank
            real(c_double), value :: eps
            integer(c_int) :: rank

        end subroutine
    end interface

    interface
        subroutine get_lr(lowrank, l, r) &
        bind(c, name = "get_lr_c")

            use iso_c_binding

            type(c_ptr) :: lowrank
            real(c_double) :: l(*)
            real(c_double) :: r(*)

        end subroutine
    end interface

    interface
        subroutine initialize_hodlr(N, M, eps, kernel, lowrank_method, is_sym, is_pd, hodlr) &
        bind(c, name = "initialize_hodlr_c")
            
            use iso_c_binding

            integer(c_int), value :: N
            integer(c_int), value :: M
            real(c_double), value :: eps
            type(c_ptr) :: kernel
            character(c_char) :: lowrank_method(*)
            logical(c_bool), value :: is_sym
            logical(c_bool), value :: is_pd
            type(c_ptr) :: hodlr

        end subroutine
    end interface

    interface
        subroutine matmat_product(hodlr, x, b) &
        bind(c, name = "matmat_product_c")

            use iso_c_binding

            type(c_ptr) :: hodlr
            real(c_double) :: x(*)
            real(c_double) :: b(*)

        end subroutine
    end interface

    interface
        subroutine factorize(hodlr) &
        bind(c, name = "factorize_c")

            use iso_c_binding
            type(c_ptr) :: hodlr

        end subroutine
    end interface

    interface
        subroutine solve(hodlr, b, x) &
        bind(c, name = "solve_c")

            use iso_c_binding

            type(c_ptr) :: hodlr
            real(c_double) :: b(*)
            real(c_double) :: x(*)

        end subroutine
    end interface

    interface
        subroutine logdeterminant(hodlr, log_det) &
        bind(c, name = "logdeterminant_c")

            use iso_c_binding

            type(c_ptr) :: hodlr
            real(c_double) :: log_det

        end subroutine
    end interface

    interface
        subroutine symmetric_factor_transpose_product(hodlr, x, b) &
        bind(c, name = "symm_factor_transpose_prod_c")

            use iso_c_binding

            type(c_ptr) ::hodlr
            real(c_double) :: x(*)
            real(c_double) :: b(*)

        end subroutine
    end interface

    interface
        subroutine symmetric_factor_product(hodlr, x, b) &
        bind(c, name = "symm_factor_prod_c")

            use iso_c_binding

            type(c_ptr) ::hodlr
            real(c_double) :: x(*)
            real(c_double) :: b(*)

        end subroutine
    end interface

    interface
        subroutine get_symmetric_factor(hodlr, W_matrix) &
        bind(c, name = "get_symm_factor_c")

            use iso_c_binding

            type(c_ptr) :: hodlr
            real(c_double) :: W_matrix(*)

        end subroutine
    end interface

end module hodlr_mod
