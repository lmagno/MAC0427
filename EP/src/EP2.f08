program EP2
    use        test, only: testf1, testg1, testh1, c1, dc1, d2c1, &
                           testf2, testg2, testh2, c2, dc2, d2c2, &
                           testf3, testg3, testh3, c3, dc3, d2c3

    use constrained, only: penalty
    use      stats2, only: f, g, h, &
                           setf, setg, seth, &
                           setc, setdc, setd2c, &
                           iteration, armijo, norm, angle, &
                           initstat, printheader, printstat, setmethod
implicit none
    double precision, allocatable :: x(:), x0(:)
    double precision              :: mu, gamma, sigma, theta, eps
    integer                       :: n, m

    gamma = 1e-4
    sigma = 1e-4
    theta = 1e-5
    eps   = 1e-5
    mu    = 10


    n = 2
    m = 1

    allocate(x(n))
    allocate(x0(n))
    x0(:) = 1
    call printheader()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! teste 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call initstat(n, m)
    call setmethod("teste 1")
    call setf(testf1)
    call setg(testg1)
    call seth(testh1)
    call setc(c1)
    call setdc(dc1)
    call setd2c(d2c1)

    mu    = 10

    call penalty(x, x0, n, m, mu, f, g, h, c1, dc1, d2c1, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
    print *, mu
    call printstat(x, mu)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! teste 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call initstat(n, m)
    call setmethod("teste 2")
    call setf(testf2)
    call setg(testg2)
    call seth(testh2)
    call setc(c2)
    call setdc(dc2)
    call setd2c(d2c2)

    mu    = 10

    call penalty(x, x0, n, m, mu, f, g, h, c2, dc2, d2c2, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
    print *, mu
    call printstat(x, mu)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! teste 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call initstat(n, m)
    call setmethod("teste 3")
    call setf(testf3)
    call setg(testg3)
    call seth(testh3)
    call setc(c3)
    call setdc(dc3)
    call setd2c(d2c3)

    mu    = 10

    call penalty(x, x0, n, m, mu, f, g, h, c3, dc3, d2c3, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
    print *, mu
    call printstat(x, mu)
    deallocate(x)
    deallocate(x0)
end program EP2
