module sol_chol
    use    utils, only: results
    use  entrada, only: le_sistema
    use cholesky, only: cholcol, cholrow
    use   trisys, only: forwcol, forwrow, backcol, backrow
    implicit none
contains
    function sol_chol_col(filename, res) result(status)
        character(len=*), intent(in)    :: filename
        type (results),   intent(inout) :: res

        integer :: status

        status = sol_chol_generic(filename, res, cholcol, forwcol, backcol)
    end function sol_chol_col

    function sol_chol_row(filename, res) result(status)
        character(len=*), intent(in)    :: filename
        type (results),   intent(inout) :: res

        integer :: status

        status = sol_chol_generic(filename, res, cholrow, forwrow, backrow)
    end function sol_chol_row

    ! Carrega A e b do sistema
    !     Ax = b
    ! com
    !     A ∈ ℝⁿˣⁿ
    !     b ∈ ℝⁿ
    ! definidos no arquivo filename, resolve
    ! o sistema por decomposição de Cholesky
    ! se possível com orientação dependente das
    ! funções passadas e guarda os tempos de
    ! execução e o erro do resultado em res.
    ! Retorna:
    !     0: caso a matriz seja definida positiva
    !        e o sistema foi resolvido com sucesso.
    !    -1: caso contrário.
    function sol_chol_generic(filename, res, chol, forw, back) result(status)
        character(len=*), intent(in)    :: filename
        type (results),   intent(inout) :: res
        integer                         :: chol, forw, back

        real, allocatable :: A(:, :)
        real, allocatable :: x(:), b(:)
        integer           :: n, i, status
        real              :: start, finish

        ! Interfaces para poder usar as funções
        interface
            function chol(n, A)
                integer, intent(in)    :: n
                real,    intent(inout) :: A(:, :)
            end function chol

            function forw(n, A, b, unit)
                integer, intent(in)    :: n
                real,    intent(in)    :: A(:, :)
                real,    intent(inout) :: b(:)
                logical, intent(in)    :: unit
            end function forw

            function back(n, A, b, trans)
                integer, intent(in)    :: n
                real,    intent(in)    :: A(:, :)
                real,    intent(inout) :: b(:)
                logical, intent(in)    :: trans
            end function back
        end interface

        ! Carrega o sistema
        call le_sistema(n, A, b, filename)

        ! Calcula o fator de Cholesky G da matriz A
        ! tal que A = GGᵀ
        call cpu_time(start)
        status = chol(n, A)
        call cpu_time(finish)
        res%tdecomp = finish - start

        ! Retorna a função caso a matriz não seja
        ! positiva definida
        if (status == -1) then
           return
        end if

        ! Agora temos o sistema
        !     G(Gᵀx) = b

        ! Gy = b
        call cpu_time(start)
        status = forw(n, A, b, unit = .false.)
        call cpu_time(finish)
        res%tforw = finish - start

        ! Retorna a função caso a matriz seja
        ! singular
        if (status == -1) then
           return
        end if

        ! Gᵀx = y
        call cpu_time(start)
        status = back(n, A, b, trans = .true.)
        call cpu_time(finish)
        res%tback = finish - start

        ! Retorna a função caso a matriz seja
        ! singular
        if (status == -1) then
           return
        end if

        ! Aloca e calcula o vetor solução esperado
        allocate(x(n))
        x = [(1 + mod(i-1, n/100), i = 1, n)]

        ! Calcula a norma do vetor diferença
        ! entre a solução esperada e a obtida,
        ! como forma de estimativa do erro
        res%erro = norm2(x - b)/sqrt(real(n))

        ! Libera a memória
        deallocate(x)
        deallocate(b)
        deallocate(A)
    end function sol_chol_generic
end module sol_chol
