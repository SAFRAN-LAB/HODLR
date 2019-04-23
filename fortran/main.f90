program main

    implicit none
    integer, parameter :: N = 5
    double precision, dimension(N, N) :: a

    a = generate_random_matrix(N)
    call print_matrix(N, a)

    contains
        function generate_random_matrix(N)

            integer, intent(in) :: N

            double precision :: rand
            integer :: i, j

            double precision, dimension(N, N) :: generate_random_matrix

            do i = 1, N
                do j = 1, N
                    call random_number(rand)
                    generate_random_matrix(i, j) = rand
                end do
            end do

        end function generate_random_matrix

        subroutine print_matrix(N, a)
            implicit none
            integer, intent(in) :: N
            double precision, intent(in) :: a(N, N)

            integer :: i, j

            do i = 1, N
                print *, (a(i, j), j = 1, N)
            end do
        end subroutine print_matrix

end program main
