program EP2
    use        test
    use constrained, only: penalty
    use      stats2, only: f, g, h, &
                           c, dc, d2c, &
                           setf, setg, seth, &
                           setc, setdc, setd2c, &
                           subproblem, iteration, armijo, norm, angle, &
                           initstat, printheader, printstat, setmethod
implicit none
    double precision, allocatable :: x(:), x0(:)
    double precision              :: mu, gamma, sigma, theta, eps
    integer                       :: n, m, i

    ! Chamadas do sistema para criar processos independentes
    ! Útil para matar uma conta específica sem afetar as outras
    interface
       function fork() bind(C)
         use iso_c_binding, only: c_int
         integer(c_int) :: fork
       end function fork

       function wait(i) bind(C)
         use iso_c_binding, only: c_int
         integer(c_int), VALUE :: i
         integer(c_int)        :: wait
       end function wait
    end interface

    gamma = 1e-4
    sigma = 1e-4
    theta = 1e-5
    eps   = 1e-6


    n = 2
    m = 1

    allocate(x(n))
    allocate(x0(n))
    x0(1) = 1
    x0(2) = 1
    call printheader()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! teste 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(fork() == 0) then
        call initstat(n, m)
        call setmethod("teste 1")
        call setf(testf1)
        call setg(testg1)
        call seth(testh1)
        call setc(c1)
        call setdc(dc1)
        call setd2c(d2c1)

        mu = 10
        open(19, file="data/teste1.dat", status="unknown", action="write")
        open(38, file="data/teste1_lagrangean.dat", status="unknown", action="write")
        call penalty(x, x0, n, m, mu, f, g, h, c, dc, d2c, gamma, sigma, theta, eps, subproblem, iteration, armijo, norm, angle)

        call printstat(x, mu)

        deallocate(x)
        deallocate(x0)

        close(19)
        close(38)
        call exit(0)
    end if
    i = wait(i)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! teste 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(fork() == 0) then
        call initstat(n, m)
        call setmethod("teste 2")
        call setf(testf2)
        call setg(testg2)
        call seth(testh2)
        call setc(c2)
        call setdc(dc2)
        call setd2c(d2c2)

        mu = 10

        open(19, file="data/teste2.dat", status="unknown", action="write")
        open(38, file="data/teste2_lagrangean.dat", status="unknown", action="write")
        call penalty(x, x0, n, m, mu, f, g, h, c, dc, d2c, gamma, sigma, theta, eps, subproblem, iteration, armijo, norm, angle)

        call printstat(x, mu)

        deallocate(x)
        deallocate(x0)

        close(19)
        close(38)
        call exit(0)
    end if
    i = wait(i)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! teste 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(fork() == 0) then
        m = 2
        call initstat(n, m)
        call setmethod("teste 3")
        call setf(testf3)
        call setg(testg3)
        call seth(testh3)
        call setc(c3)
        call setdc(dc3)
        call setd2c(d2c3)

        mu = 10

        open(19, file="data/teste3.dat", status="unknown", action="write")
        open(38, file="data/teste3_lagrangean.dat", status="unknown", action="write")
        call penalty(x, x0, n, m, mu, f, g, h, c, dc, d2c, gamma, sigma, theta, eps, subproblem, iteration, armijo, norm, angle)

        call printstat(x, mu)

        deallocate(x)
        deallocate(x0)

        close(19)
        close(38)
        call exit(0)
    end if
    i = wait(i)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! teste 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(fork() == 0) then
        m = 2
        call initstat(n, m)
        call setmethod("teste 4")
        call setf(testf4)
        call setg(testg4)
        call seth(testh4)
        call setc(c4)
        call setdc(dc4)
        call setd2c(d2c4)

        mu = 10

        open(19, file="data/teste4.dat", status="unknown", action="write")
        open(38, file="data/teste4_lagrangean.dat", status="unknown", action="write")
        call penalty(x, x0, n, m, mu, f, g, h, c, dc, d2c, gamma, sigma, theta, eps, subproblem, iteration, armijo, norm, angle)

        call printstat(x, mu)

        deallocate(x)
        deallocate(x0)

        close(19)
        close(38)
        call exit(0)
    end if
    i = wait(i)

    call sleep(1)

end program EP2
