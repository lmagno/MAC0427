module BuscaLinear
    use   trisys, only: backcol, forwcol
    use cholesky, only: cholcol
    implicit none

contains
    subroutine grad(x, x0, f, g, ls, gamma, eps)
        ! Entrada
        double precision, intent(in) :: x0(:)   ! Ponto inicial
        double precision, intent(in) :: eps     ! Tolerância
        double precision, intent(in) :: gamma   ! γ > 0, Constante para condição de Armijo

        ! Função, gradiente e busca linear
        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine g(gx, x)
                double precision, intent(out) :: gx(:)
                double precision, intent(in)  :: x(:)
            end subroutine g

            subroutine ls(alpha, x, d, f, gamma, gTd, xd)
                double precision, intent(inout) :: alpha
                double precision, intent(in)    :: x(:), d(:)
                double precision, intent(in)    :: gamma, gTd
                double precision, intent(inout) :: xd(:)

                interface
                    function f(x)
                        double precision, intent(in) :: x(:)
                        double precision             :: f
                    end function f
                end interface
            end subroutine ls
        end interface

        ! Saída
        double precision, intent(out) :: x(:)

        ! Variáveis locais
        integer                       :: n       ! Dimensão do problema
        double precision, allocatable :: d(:)    ! Direção
        double precision              :: alpha   ! Passo


        ! Variáveis para evitar chamadas de função desnecessárias
        double precision, allocatable :: gx(:)   ! g(x)
        double precision              :: gTd     ! gᵀ(x)d
        double precision, allocatable :: xd(:)

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate(xd(n))
        allocate(gx(n))

        ! Valores iniciais
        x  = x0
        call g(gx, x)

        alpha = 1.d0
        gTd   = -norm2(gx)
        do while (norm2(gx) > eps)
            ! Máxima descida
            d = -gx

            ! Aproxima chute inicial de α a partir do aceito anteriormente
            alpha = alpha*gTd
            gTd   = dot_product(gx, d)
            alpha = alpha/gTd

            ! Busca linear
            call ls(alpha, x, d, f, gamma, gTd, xd)

            ! Atualiza x com o novo valor
            x  = x + alpha*d
            call g(gx, x)
        end do

        ! Libera a memória
        deallocate(d)
        deallocate(xd)
        deallocate(gx)
    end subroutine grad

    subroutine newt(x, x0, f, g, h, ls, gamma, eps)
        ! Entrada
        double precision, intent(in) :: x0(:)   ! Ponto inicial
        double precision, intent(in) :: eps     ! Tolerância
        double precision, intent(in) :: gamma   ! γ > 0, Constante para condição de Armijo

        ! Função, gradiente, hessiana e busca linear
        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine g(gx, x)
                double precision, intent(out) :: gx(:)
                double precision, intent(in)  :: x(:)
            end subroutine g

            subroutine h(hx, x)
                double precision, intent(out) :: hx(:, :)
                double precision, intent(in)  :: x(:)
            end subroutine h

            subroutine ls(alpha, x, d, f, gamma, gTd, xd)
                double precision, intent(inout) :: alpha
                double precision, intent(in)    :: x(:), d(:)
                double precision, intent(in)    :: gamma, gTd
                double precision, intent(inout) :: xd(:)

                interface
                    function f(x)
                        double precision, intent(in) :: x(:)
                        double precision             :: f
                    end function f
                end interface
            end subroutine ls
        end interface

        ! Saída
        double precision, intent(out) :: x(:)

        ! Variáveis locais
        integer                       :: n       ! Dimensão do problema
        double precision, allocatable :: d(:)    ! Direção
        integer,          allocatable :: p(:)    ! Vetor de permutação
        double precision              :: alpha   ! Passo
        double precision              :: rho     ! Constante para shift nos autovalores
                                                 ! da hessiana

        integer                       :: i, status

        ! Variáveis para evitar chamadas de função desnecessárias
        double precision, allocatable :: gx(:)   ! ∇f(x)
        double precision, allocatable :: xd(:)
        double precision, allocatable :: hx(:,:) ! ∇²f(x)
        double precision, allocatable :: hr(:,:) ! ∇²f(x) + ρI
        double precision              :: gTd     ! ∇ᵀf(x)d

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate( p(n))
        allocate(gx(n))
        allocate(xd(n))
        allocate(hx(n, n))
        allocate(hr(n, n))

        ! Valores iniciais
        x  = x0
        call g(gx, x)

        do while (norm2(gx) > eps)
            ! Calcula fator de Cholesky de uma aproximação
            ! definida positiva da hessiana
            !     GGᵀ = ∇²f(x) + ρI
            call h(hx, x)

            ! Verifica se a própria hessiana não é definida positiva
            hr  = hx
            status = cholcol(n, hr)

            rho = 1e-3
            do while (status == -1)
                ! Não é definida positiva

                ! Faz shift nos autovalores e tenta novamente
                ! hr = ∇²f(x) + ρI
                hr = hx
                do i = 1, n
                    hr(i, i) = hr(i, i) + rho
                end do

                status = cholcol(n, hr)
                rho = 10*rho
            end do

            d = -gx
            ! Resolve Gᵀd solução da equação
            !    GGᵀd = -∇f(x)
            status = forwcol(n, hr, d, unit = .false.)
            if (status == -1) then
                print *, "A hessiana é singular!"
                call exit(1)
            end if

            ! Resolve d solução da equação
            !     Gᵀd = -G⁻¹∇f
            status = backcol(n, hr, d, trans = .true.)
            if (status == -1) then
                call exit(1)
            end if

            gTd = dot_product(gx, d)

            ! Busca linear
            alpha = 1.d0
            call ls(alpha, x, d, f, gamma, gTd, xd)

            ! Atualiza x com o novo valor
            x  = x + alpha*d
            call g(gx, x)
        end do

        ! Libera a memória
        deallocate(d)
        deallocate(p)
        deallocate(xd)
        deallocate(gx)
        deallocate(hx)
        deallocate(hr)
    end subroutine newt

    subroutine bfgs(x, x0, f, g, ls, gamma, eps)
        ! Entrada
        double precision, intent(in) :: x0(:)   ! Ponto inicial
        double precision, intent(in) :: eps     ! Tolerância
        double precision, intent(in) :: gamma   ! γ > 0, Constante para condição de Armijo

        ! Função, gradiente e busca linear
        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine g(gx, x)
                double precision, intent(out) :: gx(:)
                double precision, intent(in)  :: x(:)
            end subroutine g

            subroutine ls(alpha, x, d, f, gamma, gTd, xd)
                double precision, intent(inout) :: alpha
                double precision, intent(in)    :: x(:), d(:)
                double precision, intent(in)    :: gamma, gTd
                double precision, intent(inout) :: xd(:)

                interface
                    function f(x)
                        double precision, intent(in) :: x(:)
                        double precision             :: f
                    end function f
                end interface
            end subroutine ls
        end interface

        ! Saída
        double precision, intent(out) :: x(:)

        ! Variáveis locais
        integer                       :: n       ! Dimensão do problema
        double precision, allocatable :: d(:)    ! Direção
        double precision              :: alpha   ! Passo

        integer                       :: i, j, k

        ! Variáveis para evitar chamadas de função desnecessárias
        double precision, allocatable :: gx(:)   ! ∇f(x)
        double precision, allocatable :: xd(:)   ! x + αd
        double precision, allocatable :: h(:,:)  ! H(x)
        double precision, allocatable :: p(:)    !    x + αd  -    x
        double precision, allocatable :: q(:)    ! ∇f(x + αd) - ∇f(x)
        double precision              :: gTd     ! ∇ᵀf(x)d

        ! Constantes para atualizar a matriz pelo BFGS
        double precision              :: a
        double precision              :: pTq     ! pᵀq
        double precision              :: qTHq    ! qᵀHq
        double precision              :: pqTHij  ! (pqᵀH)ᵢⱼ
        double precision              :: HqpTij  ! (Hqpᵀ)ᵢⱼ

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate( p(n))
        allocate( q(n))
        allocate(gx(n))
        allocate(xd(n))
        allocate(h(n, n))

        ! Valores iniciais
        x  = x0
        call g(gx, x)

        ! h começa sendo a identidade
        h(:, :) = 0.d0
        do i = 1, n
            h(i, i) = 1.d0
        end do

        do while (norm2(gx) > eps)
            ! Direção de Quasi-Newton
            ! d = -H∇f
            d(:) = 0.d0
            do j = 1, n
                do i = 1, n
                    d(i) = d(i) - h(i, j)*gx(j)
                end do
            end do

            gTd = dot_product(gx, d)

            ! Busca linear
            alpha = 1.d0
            call ls(alpha, x, d, f, gamma, gTd, xd)

            q = gx

            ! Atualiza x
            x  = x + alpha*d
            call g(gx, x)

            p = alpha*d   ! p =    x + αd  -    x
            q = gx - q    ! q = ∇f(x + αd) - ∇f(x)
            pTq  = dot_product(p, q)

            ! Atualiza H
            qTHq = 0.d0
            do j = 1, n
                do i = 1, n
                    qTHq = qTHq + q(i)*h(i, j)*q(j)
                end do
            end do

            a = (1 + qTHq)/pTq

            ! Atualiza a parte triangular inferior de H (sem a diagonal principal)
            do j = 1, n
                do i = j+1, n
                    ! Explora o fato que
                    !     H = Hᵀ
                    pqTHij = 0.d0
                    HqpTij = 0.d0

                    ! Só acessa a parte triangular superior de H (diagonal principal inclusa)
                    !     i <= j
                    do k = 1, j
                        pqTHij = pqTHij + p(i)*q(k)*h(k, j)
                    end do

                    do k = j+1, n
                        pqTHij = pqTHij + p(i)*q(k)*h(j, k)
                    end do

                    do k = 1, j
                        HqpTij = HqpTij + h(i, k)*q(k)*p(j)
                    end do

                    do k = j+1, n
                        HqpTij = HqpTij + h(k, i)*q(k)*p(j)
                    end do

                    ! i > j
                    h(i, j) = h(i, j) + (a*p(i)*p(j) - pqTHij - HqpTij)/pTq
                end do
            end do

            ! Atualiza a diagonal principal
            ! Aqui é importante o fato de que o termo Hᵢᵢ só é usado
            ! na iteração i e, portanto, pode ser alterado sem afetar as seguintes
            do i = 1, n
                pqTHij = 0.d0
                HqpTij = 0.d0

                ! Só acessa a parte triangular superior de H (diagonal principal inclusa)
                !     i <= j
                do k = 1, i
                    pqTHij = pqTHij + p(i)*q(k)*h(k, i)
                end do

                do k = i+1, n
                    pqTHij = pqTHij + p(i)*q(k)*h(i, k)
                end do

                do k = 1, i
                    HqpTij = HqpTij + h(i, k)*q(k)*p(i)
                end do

                do k = i+1, n
                    HqpTij = HqpTij + h(k, i)*q(k)*p(i)
                end do

                ! i = j
                h(i, i) = h(i, i) + (a*p(i)*p(i) - pqTHij - HqpTij)/pTq
            end do

            ! Atualiza a parte triangular superior de H
            do j = 1, n
                do i = 1, j-1
                    ! i < j
                    h(i, j) = h(j, i)
                end do
            end do
        end do

        ! Libera a memória
        deallocate(d)
        deallocate(p)
        deallocate(q)
        deallocate(xd)
        deallocate(h)
        deallocate(gx)
    end subroutine bfgs
end module BuscaLinear
