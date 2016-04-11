module cholesky
    implicit none
contains
    ! Calcula o fator de Cholesky pelo método do produto externo,
    ! orientado a colunas.
    ! Guarda o fator na parte triangular inferior da matriz A e
    ! retorna:
    !     0: se o fator foi calculado com sucesso
    !    -1: se a matriz A não é positiva definida e não foi
    !        possível calcular seu fator
    function cholcol(n, A) result(status)
        integer, intent(in)    :: n
        real,    intent(inout) :: A(:, :)

        real    :: akk
        integer :: status
        integer :: i, j, k

        status = -1
        do k = 1, n
            akk = A(k, k)
            if (akk <= 0) then
                return
            end if

            akk = sqrt(akk)
            A(k, k) = akk
            do i = k+1, n
                A(i, k) = A(i, k)/akk
            end do

            do j = k+1, n
                do i = j, n
                    A(i, j) = A(i, j) - A(i, k)*A(j, k)
                end do
            end do
        end do

        status = 0
    end function cholcol

    ! Calcula o fator de Cholesky pelo método do produto externo,
    ! orientado a linhas.
    ! Guarda o fator na parte triangular inferior da matriz A e
    ! retorna:
    !     0: se o fator foi calculado com sucesso
    !    -1: se a matriz A não é positiva definida e não foi
    !        possível calcular seu fator
    function cholrow(n, A) result(status)
        integer, intent(in)    :: n
        real,    intent(inout) :: A(:, :)

        real    :: akk
        integer :: status
        integer :: i, j, k

        status = -1
        do k = 1, n
            akk = A(k, k)
            if (akk <= 0) then
                return
            end if

            akk = sqrt(akk)
            A(k, k) = akk
            do i = k+1, n
                A(i, k) = A(i, k)/akk
            end do

            do i = k+1, n
                do j = k+1, i
                    A(i, j) = A(i, j) - A(i, k)*A(j, k)
                end do
            end do
        end do

        status = 0
    end function cholrow

end module cholesky
