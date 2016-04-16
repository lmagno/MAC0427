module LS
    implicit none

contains
    subroutine lsquad(alpha, x, d, f, gamma, gTd, xd, armijo)
        double precision, intent(inout) :: alpha
        double precision, intent(in)    :: x(:), d(:)
        double precision, intent(in)    :: gamma, gTd
        double precision, intent(inout) :: xd(:) ! Para economizar alocações

        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine armijo()
            end subroutine armijo
        end interface

        ! Variáveis locais
        double precision :: fx, fx0, fx1
        double precision :: a0, a1
        double precision :: a, b, min

        fx = f(x)

        ! Tenta o passo inicial
        a0 = alpha

        xd  = x + a0*d
        fx0 = f(xd)
        if (fx0 <= fx + a0*gamma*gTd) then
            call armijo()
            ! Condição de Armijo
            return
        end if

        a1  = a0
        fx1 = fx0
        do while (fx1 > fx + a1*gamma*gTd)
            call armijo()
            ! Aproximação quadrática de ϕ(α) = f(x + αd)
            !     q(α) = aα² + bα + q(0)
            ! cujo minimizador é
            !     α = -b/2a

            a = (fx0 - fx - a0*gTd)/(a0*a0)
            b = gTd

            min = -b/(2.d0*a)

            ! Escolhe α₁ tal que 0.1α₀ <= α₁ <= 0.9α₀
            if (min < 0.1d0*a0) then
                a1 = 0.1d0*a0
            else if (min > 0.9*a0) then
                a1 = 0.9d0*a0
            else
                a1 = min
            end if

            xd  = x + a1*d
            fx1 = f(xd)

            a0  = a1
            fx0 = fx1
        end do

        alpha = a1
    end subroutine lsquad

    subroutine lscube(alpha, x, d, f, gamma, gTd, xd, armijo)
        double precision, intent(inout) :: alpha
        double precision, intent(in)    :: x(:), d(:)
        double precision, intent(in)    :: gamma, gTd
        double precision, intent(inout) :: xd(:) ! Para economizar alocações

        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine armijo()
            end subroutine armijo
        end interface

        ! Variáveis locais
        double precision :: fx, fx0, fx1, fx2
        double precision :: a0, a1, a2
        double precision :: a, b, min

        fx = f(x)

        ! Tenta o passo inicial
        a0 = alpha

        xd  = x + a0*d
        fx0 = f(xd)
        if (fx0 <= fx + a0*gamma*gTd) then
            call armijo()
            ! Condição de Armijo
            return
        end if

        ! Aproximação quadrática de ϕ(α) = f(x + αd)
        !     q(α) = aα² + bα + q(0)
        ! cujo minimizador é
        !     α = -b/2a

        a = (fx0 - fx - a0*gTd)/(a0*a0)
        b = gTd

        min = -b/(2.d0*a)

        ! Escolhe α₁ tal que 0.1α₀ <= α₁ <= 0.9α₀
        if (min < 0.1d0*a0) then
            a1 = 0.1d0*a0
        else if (min > 0.9*a0) then
            a1 = 0.9d0*a0
        else
            a1 = min
        end if

        xd  = x + a1*d
        fx1 = f(xd)
        if (fx1 <= fx + a1*gamma*gTd) then
            call armijo()
            alpha = a1
            return
        end if

        a2  = a1
        fx2 = fx1
        do while (fx2 > fx + a2*gamma*gTd)
            call armijo()
            ! Aproximação cúbica de ϕ(α) = f(x + αd)
            !     q(α) = aα³ + bα² + q'(0)α + q(0)
            ! cujo minimizador é
            !     α = (-b + √(b² - 3aq'(0)))/3a
            ! onde
            !     q(0) = fx
            !    q'(0) = gTd

            a = (    a0*a0*(fx1 - fx - gTd*a1) -    a1*a1*(fx0 - fx - gTd*a0))/(a0*a0*a1*a1*(a1-a0))
            b = (-a0*a0*a0*(fx1 - fx - gTd*a1) + a1*a1*a1*(fx0 - fx - gTd*a0))/(a0*a0*a1*a1*(a1-a0))

            min = (-b + sqrt(b*b - 3*a*gTd))/(3*a)

            ! Escolhe α₂ tal que 0.1α₁ <= α₂ <= 0.9α₁
            if (min < 0.1d0*a1) then
                a2 = 0.1d0*a1
            else if (min > 0.9*a1) then
                a2 = 0.9d0*a1
            else
                a2 = min
            end if

            xd  = x + a2*d
            fx2 = f(xd)

            a0  = a1
            a1  = a2
            fx0 = fx1
            fx1 = fx2
        end do

        alpha = a2
    end subroutine lscube
end module LS
