program EP1
    use     methods, only: grad, newt, bfgs
    use         mgh, only: setprob, getdim, gettries, getname, getinit, mgh_f, mgh_g, mgh_h
    use          ls, only: lsquad, lscube
    use  rosenbrock, only: rosef, roseg, roseh
    use  paraboloid, only: paraf, parag, parah
    use       stats, only: f, g, h, &
                           setf, setg, seth, &
                           iteration, armijo, norm, angle, &
                           initstat, printheader, printstat, setmethod
    implicit none

    integer                       :: p, s
    integer                       :: n, nprob, ntry, ntries, i
    character(len=30)             :: name
    double precision              :: factor = 1.d0
    double precision, allocatable :: x(:), x0(:)
    double precision              :: gamma, sigma, theta, eps

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
    eps   = 1e-5

    call printheader()

    !!!!!!!!!!!!!!!!!!!!!!!!!!! Parabolóide !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, "Paraboloid"
    call initstat()
    call setf(paraf)
    call setg(parag)
    call seth(parah)

    n = 2
    allocate(x(2))
    allocate(x0(2))
    x0(1) = 100
    x0(2) = 100

    INCLUDE 'src/runmethods.f08'
    deallocate(x)
    deallocate(x0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!! Rosenbrock !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, "Rosenbrock"
    call initstat()
    call setf(rosef)
    call setg(roseg)
    call seth(roseh)

    n = 2
    allocate(x(2))
    allocate(x0(2))
    x0(1) = 1
    x0(2) = 2

    INCLUDE 'src/runmethods.f08'
    deallocate(x)
    deallocate(x0)

    !!!!!!!!!!!!!!!!!!!!!!! Moré-Garbow-Hillstrom !!!!!!!!!!!!!!!!!!!!!!!!!
    call initstat()
    call setf(mgh_f)
    call setg(mgh_g)
    call seth(mgh_h)

    do nprob = 1, 18
        call setprob(nprob)
        ntries = gettries()
        factor = 1.d0

        name = getname()
        n    = getdim()

        print '(i2, x, a)', nprob, name
        allocate( x(n))
        allocate(x0(n))
        call getinit(x0, factor)

        INCLUDE 'src/runmethods.f08'

        deallocate(x)
        deallocate(x0)
     end do
end program EP1
