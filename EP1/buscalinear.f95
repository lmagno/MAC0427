module BuscaLinear
    implicit none

contains
    subroutine maximaDescida(x, x0, f, g, gamma, eps)
        ! Entrada
        real, intent(in) :: x0(:)   ! Ponto inicial
        real, intent(in) :: eps     ! Tolerância
        real, intent(in) :: gamma   ! γ > 0, Constante para condição de Armijo

        ! Função e suas derivadas
        interface
            function f(x)
                real, intent(in) :: x(:)
                real             :: f
            end function f

            function g(x)
                real, intent(in) :: x(:)
                real             :: g(size(x))
            end function g
        end interface

        ! Saída
        real, intent(out) :: x(:)

        ! Variáveis locais
        integer           :: n       ! Dimensão do problema
        real, allocatable :: d(:)    ! Direção
        real              :: alpha   ! Passo
        real              :: a, b    ! Coeficientes da aproximação quadrática
                                     ! de f

        ! Variáveis para evitar chamadas de função desnecessárias
        real              :: fx      ! f(x)
        real, allocatable :: gx(:)   ! g(x)
        real              :: gTd     ! gᵀ(x)d
        real              :: fxd     ! f(x + alpha*d)

        ! Aloca
        n = size(x0)
        allocate( d(n))
        allocate(gx(n))

        ! Valores iniciais
        x  = x0
        fx = f(x)
        gx = g(x)

        do while (norm2(gx) > eps)
            print *, norm2(gx)
            print *, x
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Máxima descida
            d = -gx
            gTd = dot_product(gx, d)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Passo 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            alpha = 1.0

            ! Condição de Armijo
            fxd = f(x + alpha*d)
            do while (fxd > fx + alpha*gamma*gTd)
                print *, fxd - fx - alpha*gamma*gTd
                ! Aproximação quadrática de f perto de x ao longo de d:
                !     q(alpha) = a*alpha^2 + b*alpha + c
                ! cujo minimizador é
                !     alpha = -b/2a

                a = (fxd - fx - alpha*gTd)/(alpha*alpha)
                b = gTd

                alpha = -b/(2*a)
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
end module BuscaLinear
