module BuscaLinear
    use  utils, only: Results
    use     lu, only: lucol, sscol
    implicit none

contains
    subroutine grad(x, x0, f, g, ls, gamma, eps)
        ! Entrada
        double precision, intent(in) :: x0(:)   ! Ponto inicial
        double precision, intent(in) :: eps     ! Tolerância
        double precision, intent(in) :: gamma   ! γ > 0, Constante para condição de Armijo

        ! Função e suas derivadas
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
        double precision              :: a, b    ! Coeficientes da aproximação quadrática
                                                 ! de f
        double precision :: min

        ! Variáveis para evitar chamadas de função desnecessárias
        double precision              :: fx      ! f(x)
        double precision, allocatable :: gx(:)   ! g(x)
        double precision              :: gTd     ! gᵀ(x)d
        double precision, allocatable :: xd(:)
        double precision              :: fxd     ! f(x + alpha*d)

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate(xd(n))
        allocate(gx(n))

        ! Valores iniciais
        x  = x0
        fx = f(x)
        call g(gx, x)

        alpha = 1.d0
        gTd   = -norm2(gx)
        do while (norm2(gx) > eps)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            fx = f(x)
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

        ! Função e suas derivadas
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
        double precision              :: a, b    ! Coeficientes da aproximação quadrática
                                                 ! de f
        double precision              :: min
        type(Results)                 :: res
        integer                       :: status

        ! Variáveis para evitar chamadas de função desnecessárias
        double precision              :: fx      ! f(x)
        double precision, allocatable :: gx(:)   ! g(x)
        double precision, allocatable :: xd(:)
        double precision, allocatable :: hx(:,:) ! h(x)
        double precision              :: gTd     ! gᵀ(x)d
        double precision              :: fxd     ! f(x + alpha*d)

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate( p(n))
        allocate(gx(n))
        allocate(xd(n))
        allocate(hx(n, n))

        ! Valores iniciais
        x  = x0
        fx = f(x)
        call g(gx, x)

        do while (norm2(gx) > eps)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Direção de Newton
            call h(hx, x)
            d = -gx
            status = lucol(n, hx, p)
            if (status == -1) then
                print *, "A hessiana é singular!"
                call exit(1)
            end if

            status = sscol(n, hx, p, d, res)
            if (status == -1) then
                call exit(1)
            end if

            gTd = dot_product(gx, d)

            ! Busca linear
            alpha = 1.d0
            call ls(alpha, x, d, f, gamma, gTd, xd)

            ! Atualiza x com o novo valor
            x  = x + alpha*d
            fx = f(x)
            call g(gx, x)
        end do

        ! Libera a memória
        deallocate(d)
        deallocate(p)
        deallocate(xd)
        deallocate(hx)
        deallocate(gx)
    end subroutine newt

    subroutine bfgs(x, x0, f, g, gamma, eps)
        ! Entrada
        double precision, intent(in) :: x0(:)   ! Ponto inicial
        double precision, intent(in) :: eps     ! Tolerância
        double precision, intent(in) :: gamma   ! γ > 0, Constante para condição de Armijo

        ! Função e suas derivadas
        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine g(gx, x)
                double precision, intent(out) :: gx(:)
                double precision, intent(in)  :: x(:)
            end subroutine g
        end interface

        ! Saída
        double precision, intent(out) :: x(:)

        ! Variáveis locais
        integer                       :: n       ! Dimensão do problema
        double precision, allocatable :: d(:)    ! Direção
        double precision              :: alpha   ! Passo
        double precision              :: a, b    ! Coeficientes da aproximação quadrática
                                                 ! de f
        double precision              :: min
        type(Results)                 :: res
        integer                       :: i, j, k

        ! Variáveis para evitar chamadas de função desnecessárias
        double precision              :: fx      ! f(x)
        double precision, allocatable :: gx(:)   ! g(x)
        double precision, allocatable :: xd(:)   ! x + αd
        double precision              :: fxd     ! f(x + αd)
        double precision, allocatable :: h(:,:)  ! h(x)
        double precision, allocatable :: p(:)    ! x + αd    - x
        double precision, allocatable :: q(:)    ! g(x + αd) - g(x)
        double precision              :: gTd     ! gᵀ(x)d

        ! Constantes para atualizar a matriz pelo BFGS
        double precision              :: pTq     ! pᵀq
        double precision              :: qTHq    ! qᵀHq
        double precision              :: pqTHij  ! (pqᵀH)ᵢⱼ
        double precision              :: HqpTij  ! (Hqpᵀ)ᵢⱼ

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate( p(n))
        allocate(gx(n))
        allocate(xd(n))
        allocate(h(n, n))

        ! Valores iniciais
        x  = x0
        fx = f(x)
        call g(gx, x)

        ! hx começa sendo a identidade
        h(:, :) = 0.d0
        do i = 1, n
            h(i, i) = 1.d0
        end do

        d(:) = 0.d0
        do while (norm2(gx) > eps)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Direção de Quasi-Newton
            ! d = -Hg(x)
            do j = 1, n
                do i = 1, n
                    d(i) = d(i) - h(i, j)*gx(j)
                end do
            end do
            gTd = dot_product(gx, d)

            ! Condição de Armijo
            alpha = 1.d0

            xd  = x + alpha*d
            fxd = f(xd)
            do while (fxd > fx + alpha*gamma*gTd)
                ! Aproximação quadrática de f perto de x ao longo de d:
                !     q(alpha) = a*alpha^2 + b*alpha + c
                ! cujo minimizador é
                !     alpha = -b/2a

                a = (fxd - fx - alpha*gTd)/(alpha*alpha)
                b = gTd

                min = -b/(2.d0*a)

                ! Escolhe novo α tal que 0.1α <= α <= 0.9α
                if (min < 0.1d0*alpha) then
                    alpha = 0.1d0*alpha
                else if (min > 0.9*alpha) then
                    alpha = 0.9d0*alpha
                else
                    alpha = min
                end if

                xd  = x + alpha*d
                fxd = f(xd)
            end do
            ! print *, alpha, norm2(d)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            p  = alpha*d
            q  = gx

            ! Atualiza x
            x  = xd
            fx = fxd
            call g(gx, x)

            ! q = g(x+αd) - g(x)
            q  = gx - q
            ! print *, norm2(gx), norm2(p), norm2(q)

            pTq  = dot_product(p, q)
            qTHq = 0.d0
            do j = 1, n
                do i = 1, n
                    qTHq = qTHq + q(i)*h(i, j)*q(j)
                end do
            end do

            a = (1.d0 + qTHq)/pTq
            do j = 1, n
                do i = 1, n
                    pqTHij = 0.d0
                    HqpTij = 0.d0

                    do k = 1, n
                        pqTHij = pqTHij + p(i)*q(k)*h(k, j)
                    end do

                    do k = 1, n
                        HqpTij = HqpTij + h(i, k)*q(k)*p(j)
                    end do

                    h(i, j) = h(i, j) + (a*p(i)*p(j) - pqTHij - HqpTij)/pTq
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
