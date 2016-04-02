program EP1
    use BuscaLinear, only: maximaDescida
    use MGH,         only: setprob, getdim, getname, getinit, f, g
    implicit none

    integer                       :: n
    character(len=30)             :: name
    double precision              :: factor = 1.d0
    double precision, allocatable :: x(:), x0(:)
    double precision              :: gamma, eps

    call setprob(1)
    name = getname()
    n    = getdim()
    x0   = getinit(factor)
    allocate(x(n))

    print *, name
    print *, "n     =", n
    print *, "x0    =", x
    print *, "f(x0) =", f(x)
    print *, "g(x0) =", g(x)

    gamma = 0.5
    eps   = 1e-6

    call maximaDescida(x, x0, f, g, gamma, eps)

    print *, ""
    print *, "Mínimo:"
    print *, "x    =", x
    print *, "f(x) =", f(x)
    print *, "g(x) =", g(x)

    deallocate(x)
    deallocate(x0)
end program EP1
