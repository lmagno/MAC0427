module rosenbrock
    implicit none

contains
    function rosef(x)
        double precision, intent(in) :: x(:)
        double precision             :: rosef

        double precision :: x1, x2
        double precision :: a, b

        x1 = x(1)
        x2 = x(2)

        a = (1 - x1)
        b = (x2 - x1*x1)

        rosef = a*a + 100*b*b
    end function rosef

    subroutine roseg(gx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: gx(:)

        double precision :: x1, x2

        x1 = x(1)
        x2 = x(2)

        gx(1) = 400*x1*x1*x1 + 2*x1 -2 -400*x1*x2
        gx(2) = -200*x1*x1 + 200*x2
    end subroutine roseg

    subroutine roseh(hx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: hx(:, :)

        double precision :: x1, x2

        x1 = x(1)
        x2 = x(2)

        hx(1, 1) = 1200*x1*x1 + 2 - 400*x2
        hx(1, 2) = -400*x1
        hx(2, 1) = -400*x1
        hx(2, 2) = 200
    end subroutine roseh
end module rosenbrock

program EP1
    use BuscaLinear, only: grad, newt, bfgs
    use         MGH, only: setprob, getdim, gettries, getname, getinit, mgh_f, mgh_g, mgh_h
    use          ls, only: lsquad, lscube
    use  rosenbrock, only: rosef, roseg, roseh
    use       stats, only: f, g, h, &
                           setf, setg, seth, &
                           iteration, armijo, norm, angle, &
                           initstat, printheader, printstat, setmethod, sett0, settf
    implicit none

    integer                       :: p, s
    integer                       :: n, nprob, ntry, ntries, i
    character(len=30)             :: name
    double precision              :: factor = 1.d0
    double precision, allocatable :: x(:), x0(:)
    double precision              :: gamma, sigma, theta, eps

    double precision              :: tf, t0
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


    call initstat()
    call setf(mgh_f)
    call setg(mgh_g)
    call seth(mgh_h)

    call printheader()
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

        ! call initstat()
        ! call setf(rosef)
        ! call setg(roseg)
        ! call seth(roseh)
        !
        ! n = 2
        ! allocate(x(2))
        ! allocate(x0(2))
        ! x0(1) = 1
        ! x0(2) = 2

        ! Cria novo processo
        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("gradquad")

            call sett0()
            call grad(x, x0, f, g, lsquad, gamma, eps, iteration, armijo)
            call settf()

            call printstat(x)

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("gradcube")

            call sett0()
            call grad(x, x0, f, g, lscube, gamma, eps, iteration, armijo)
            call settf()

            call printstat(x)

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("newtquad")

            call sett0()
            call newt(x, x0, f, g, h, lsquad, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
            call settf()

            call printstat(x)


            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("newtcube")

            call sett0()
            call newt(x, x0, f, g, h, lscube, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
            call settf()

            call printstat(x)

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("bfgsquad")

            call sett0()
            call bfgs(x, x0, f, g, lsquad, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
            call settf()

            call printstat(x)


            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("bfgscube")

            call sett0()
            call bfgs(x, x0, f, g, lscube, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
            call settf()

            call printstat(x)

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        ! Espera os processos filhos terminarem
        do i = 1, 6
            p = wait(s)
        end do

        deallocate(x)
        deallocate(x0)
     end do
end program EP1
