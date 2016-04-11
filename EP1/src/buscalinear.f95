module BuscaLinear
    use  utils, only: Results
    use     lu, only: lucol, sscol
    implicit none

contains
    subroutine maximaDescida(x, x0, f, g, gamma, eps)
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

            function g(x)
                double precision, intent(in) :: x(:)
                double precision             :: g(size(x))
            end function g
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
        double precision              :: fxd     ! f(x + alpha*d)

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate(gx(n))

        ! Valores iniciais
        x  = x0
        fx = f(x)
        gx = g(x)

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

            ! Condição de Armijo
            fxd = f(x + alpha*d)
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

                fxd = f(x + alpha*d)
            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            x  = x + alpha*d
            fx = f(x)
            gx = g(x)
        end do

        ! Libera a memória
        deallocate(d)
        deallocate(gx)
    end subroutine maximaDescida

    subroutine newton(x, x0, f, g, h, gamma, eps)
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

            function g(x)
                double precision, intent(in) :: x(:)
                double precision             :: g(size(x))
            end function g

            function h(x)
                double precision, intent(in) :: x(:)
                double precision             :: h(size(x), size(x))
            end function h
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
        double precision, allocatable :: hx(:,:) ! h(x)
        double precision              :: gTd     ! gᵀ(x)d
        double precision              :: fxd     ! f(x + alpha*d)

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate( p(n))
        allocate(gx(n))
        allocate(hx(n, n))

        ! Valores iniciais
        x  = x0
        fx = f(x)
        gx = g(x)

        do while (norm2(gx) > eps)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Direção de Newton
            hx = h(x)
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

            ! Condição de Armijo
            alpha = 1.d0
            fxd = f(x + alpha*d)
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

                fxd = f(x + alpha*d)
            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            x  = x + alpha*d
            fx = f(x)
            gx = g(x)
        end do

        ! Libera a memória
        deallocate(d)
        deallocate(p)
        deallocate(hx)
        deallocate(gx)
    end subroutine newton

end module BuscaLinear
