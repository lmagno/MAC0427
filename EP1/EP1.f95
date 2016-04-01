module Foo
    implicit none
contains

    function f(x)
        real, intent(in) :: x(:)
        real             :: f

        f = (10*(x(2) - x(1)*x(1)))**2 + (1 - x(1))**2
    end function f

    function g(x)
        real, intent(in) :: x(:)
        real             :: g(size(x))

        g(1) = 2*(10*(x(2) - x(1)*x(1)))*(-20*x(1)) + 2*(1 - x(1))*(-1)
        g(2) = 2*(10*(x(2) - x(1)*x(1)))*(10)
    end function g

end module Foo

program EP1
    use BuscaLinear, only: maximaDescida
    use Foo,         only: f, g
    implicit none

    integer           :: n
    real, allocatable :: x(:), x0(:)
    real              :: gamma, eps

    n = 2
    allocate( x(n))
    allocate(x0(n))

    x0(1) = -1.2
    x0(2) = 1
    gamma = 0.5
    eps   = 1e-6

    call maximaDescida(x, x0, f, g, gamma, eps)

    print *, x
    deallocate(x)
    deallocate(x0)
end program EP1
