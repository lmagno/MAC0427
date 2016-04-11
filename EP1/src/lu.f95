module lu
    use  utils, only: swap, Results
    use trisys, only: forwcol, forwrow, backcol, backrow
    implicit none
contains
    ! Dada uma matriz A ∈ ℝⁿˣⁿ, calcula
    ! sua fatoração LU tal que
    !     PA = LU
    ! onde
    !     P é uma matriz de permutação
    !     L é uma matriz triangular inferior unitária
    !     U é uma matriz triangular superior
    ! com orientação a colunas.
    !    Caso A não seja singular, guarda
    ! L na parte triangular inferior (sem a diagonal
    ! principal) de A e U na parte triangular superior.
    !    Além disso, guarda no vetor p uma
    ! representação da matriz P.
    ! Retorna:
    !     0: se a matriz não for singular e a fatoração
    !        foi completada com sucesso
    !    -1: se a matriz for singular
    function lucol(n, A, p) result(status)
        integer,          intent(in)    :: n
        double precision, intent(inout) :: A(:, :)
        integer,          intent(out)   :: p(:)

        integer :: i, j, k
        integer :: imax, status

        status = -1
        do k = 1, n

            ! Acha o máximo elemento absoluto da coluna k
            ! e guarda a linha em que ele aparece
            imax = k
            do i = k+1, n
                if (abs(A(i, k)) > abs(A(imax, k))) then
                    imax = i
                end if
            end do

            ! Verifica se a matriz é singular
            if (A(imax, k) == 0.0) then
                return
            end if

            p(k) = imax

            ! Troca as linhas de forma a deixar o maior elemento
            ! da coluna k na posição de pivô
            if (imax /= k) then
                do j = 1, n
                    call swap(A(k, j), A(imax, j))
                end do
            end if

            ! Calcula os multiplicadores
            do i = k+1, n
                A(i, k) = A(i, k)/A(k, k)
            end do

            ! Reduz as linhas
            do j = k+1, n
                do i = k+1, n
                    A(i, j) = A(i, j) - A(k, j)*A(i, k)
                end do
            end do
        end do

        status = 0
    end function lucol

    ! Dada uma matriz A ∈ ℝⁿˣⁿ, calcula
    ! sua fatoração LU tal que
    !     PA = LU
    ! onde
    !     P ∈ ℝⁿˣⁿ é uma matriz de permutação
    !     L ∈ ℝⁿˣⁿ é uma matriz triangular inferior unitária
    !     U ∈ ℝⁿˣⁿ é uma matriz triangular superior
    ! com orientação a colunas.
    !    Caso A não seja singular, guarda
    ! L na parte triangular inferior (sem a diagonal
    ! principal) de A e U na parte triangular superior.
    !    Além disso, guarda no vetor p uma
    ! representação da matriz P.
    ! Retorna:
    !     0: se a matriz não for singular e a fatoração
    !        foi completada com sucesso
    !    -1: se a matriz for singular
    function lurow(n, A, p) result(status)
        integer, intent(in)    :: n
        double precision,    intent(inout) :: A(:, :)
        integer, intent(out)   :: p(:)

        integer :: i, j, k
        integer :: imax, status

        status = -1
        do k = 1, n

            ! Acha o máximo elemento absoluto da coluna k
            ! e guarda a linha em que ele aparece
            imax = k
            do i = k+1, n
                if (abs(A(i, k)) > abs(A(imax, k))) then
                    imax = i
                end if
            end do

            ! Verifica se a matriz é singular
            if (A(imax, k) == 0.0) then
                return
            end if

            p(k) = imax

            ! Troca as linhas de forma a deixar o maior elemento
            ! da coluna k na posição de pivô
            if (imax /= k) then
                do j = 1, n
                    call swap(A(k, j), A(imax, j))
                end do
            end if

            ! Calcula os multiplicadores e
            ! reduz as linhas
            do i = k+1, n
                A(i, k) = A(i, k)/A(k, k)
                do j = k+1, n
                    A(i, j) = A(i, j) - A(k, j)*A(i, k)
                end do
            end do
        end do

        status = 0
    end function lurow

    function sscol(n, A, p, b, res) result(status)
        integer,          intent(in)    :: n, p(:)
        double precision, intent(in)    :: A(:, :)
        double precision, intent(inout) :: b(:)
        type (Results),   intent(inout) :: res

        integer :: status

        status = ss_generic(n, A, p, b, res, forwcol, backcol)
    end function sscol

    function ssrow(n, A, p, b, res) result(status)
        integer,          intent(in)    :: n, p(:)
        double precision, intent(in)    :: A(:, :)
        double precision, intent(inout) :: b(:)
        type (Results),   intent(inout) :: res

        integer :: status

        status = ss_generic(n, A, p, b, res, forwrow, backrow)
    end function ssrow

    ! Resolve o sistema (encontra x)
    !    LUx = Pb
    ! com
    !    L ∈ ℝⁿˣⁿ triangular inferior unitária
    !    U ∈ ℝⁿˣⁿ triangular superior
    !    P ∈ ℝⁿˣⁿ matriz de permutação
    !    x, b ∈ ℝⁿ
    ! com orientação dependente das funções forw
    ! e back passadas.
    !    L e U são ambas guardadas em A, e no vetor p
    ! é guardada uma representação da matriz P.
    !    Se L e U não forem singulares, calcula x
    ! e o guarda no vetor b, guardando também os tempos
    ! de execução dos passos em res.
    ! Retorna:
    !     0: caso L e U não forem singulares e o sistema
    !        seja resolvido com sucesso.
    !    -1: caso contrário.
    function ss_generic(n, A, p, b, res, forw, back) result(status)
        integer,          intent(in)    :: n, p(:)
        double precision, intent(in)    :: A(:, :)
        double precision, intent(inout) :: b(:)
        type (Results),   intent(inout) :: res
        integer                         :: forw, back

        integer          :: i, status
        double precision :: start, finish

        ! Interfaces para poder usar as funções
        interface
            function forw(n, A, b, unit)
                integer,          intent(in)    :: n
                double precision, intent(in)    :: A(:, :)
                double precision, intent(inout) :: b(:)
                logical,          intent(in)    :: unit
            end function forw

            function back(n, A, b, trans)
                integer,          intent(in)    :: n
                double precision, intent(in)    :: A(:, :)
                double precision, intent(inout) :: b(:)
                logical,          intent(in)    :: trans
            end function back
        end interface

        ! Calcula Pb
        do i = 1, n
            call swap(b(i), b(p(i)))
        end do

        ! Agora temos o sistema
        !    LUx = Pb

        ! Primeiro resolvemos
        !    Ly = Pb

        call cpu_time(start)
        status = forw(n, A, b, unit = .true.)
        call cpu_time(finish)
        res%tforw = finish - start

        if (status == -1) then
            return
        end if

        ! Agora resolvemos
        !    Ux = y
        call cpu_time(start)
        status = back(n, A, b, trans = .false.)
        call cpu_time(finish)
        res%tback = finish - start

        if (status == -1) then
            return
        end if
    end function ss_generic
end module lu
