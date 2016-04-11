! Módulo com algumas utilidades para matrizes/vetores
module utils
    implicit none

    type Results
        sequence
        double precision :: erro
        double precision :: tdecomp, tforw, tback
    end type Results

contains
    ! Imprime um vetor em formato de coluna
    subroutine pvec(v)
        double precision, intent(in) :: v(:)

        integer :: n, i

        n = size(v)

        do i = 1, n
            print *, v(i)
        end do
    end subroutine pvec

    ! Imprime uma matriz
    subroutine pmatriz(A)
        double precision, intent(in) :: A(:, :)

        integer :: s(2), i

        s = shape(A)

        do i = 1, s(1)
            print *, A(i, :)
        end do
    end subroutine pmatriz

    subroutine swap(a, b)
        double precision, intent(inout) :: a, b

        double precision :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine swap

end module utils
