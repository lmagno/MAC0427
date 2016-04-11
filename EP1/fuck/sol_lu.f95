module sol_lu
    use   utils, only: pmatriz, pvec, swap, Results
    use      lu, only: lucol, lurow, sscol, ssrow
    use entrada, only: le_sistema
    implicit none
contains
    function sol_lu_col(filename, res) result(status)
        character(len=*), intent(in)    :: filename
        type (Results),   intent(inout) :: res

        integer :: status

        status = sol_lu_generic(filename, lucol, sscol, res)
    end function sol_lu_col

    function sol_lu_row(filename, res) result(status)
        character(len=*), intent(in)    :: filename
        type (Results),   intent(inout) :: res

        integer :: status

        status = sol_lu_generic(filename, lurow, ssrow, res)
    end function sol_lu_row

    ! Carrega A e b do sistema
    !     Ax = b
    ! com
    !     A ∈ ℝⁿˣⁿ
    !     b ∈ ℝⁿ
    ! definidos no arquivo filename, resolve
    ! o sistema por decomposição LU se possível
    ! com orientação dependente das funções passadas.
    ! e guarda os tempos de execução e o erro do
    ! resultado em res.
    ! Retorna:
    !     0: caso a matriz não seja singular
    !        e o sistema foi resolvido com sucesso.
    !    -1: caso contrário.
    function sol_lu_generic(filename, lu, ss, res) result(status)
        character(len=*), intent(in)    :: filename
        type (Results),   intent(inout) :: res
        integer                         :: lu, ss

        real,    allocatable :: A(:, :)
        real,    allocatable :: x(:), b(:)
        integer, allocatable :: p(:)
        integer              :: i, n, status
        real                 :: start, finish

        ! Interfaces para poder usar as funções
        interface
           function lu(n, A, p)
             integer, intent(in)    :: n
                real, intent(inout) :: A(:, :)
             integer, intent(out)   :: p(:)
           end function lu

           function ss(n, A, p, b, res)
             use utils, only: Results
                    integer, intent(in)    :: n, p(:)
                       real, intent(in)    :: A(:, :)
                       real, intent(inout) :: b(:)
             type (Results), intent(inout) :: res
           end function ss
        end interface

        call le_sistema(n, A, b, filename)
        allocate(p(n))

        call cpu_time(start)
        status = lu(n, A, p)
        call cpu_time(finish)
        res%tdecomp = finish - start
        if (status == -1) then
           return
        end if

        status = ss(n, A, p, b, res)
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

        deallocate(A)
        deallocate(b)
        deallocate(x)
        deallocate(p)
    end function sol_lu_generic
end module sol_lu
