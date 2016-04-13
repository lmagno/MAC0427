program EP1
    use BuscaLinear, only: grad, newt, bfgs
    use         MGH, only: setprob, setmethod, getdim, gettries, getname, getinit, f, g, h, nfev, ngev, nhev
    use          ls, only: lsquad, lscube
    implicit none

    integer                       :: p, s
    integer                       :: n, nprob, ntry, ntries, i
    character(len=30)             :: name
    double precision              :: factor = 1.d0
    double precision, allocatable :: x(:), x0(:)
    double precision              :: gamma, eps

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

    print '(18x, a16, a11, a11, a13, a14)', "‖∇f(x)‖", "t (s)", "nº f", "nº ∇f", "nº ∇²f"
    gamma = 1e-4
    ! eps   = epsilon(0.d0)
    eps = 1e-7
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

        ! Cria novo processo
        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("gradquad")

            call cpu_time(t0)
            call grad(x, x0, f, g, lsquad, gamma, eps)
            call cpu_time(tf)

            call g(x0, x)
            print '(a18, e10.2, f11.3, 2i10)', "gradquad", norm2(x0), tf-t0, nfev(), ngev()

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("gradcube")

            call cpu_time(t0)
            call grad(x, x0, f, g, lscube, gamma, eps)
            call cpu_time(tf)

            call g(x0, x)
            print '(a18, e10.2, f11.3, 2i10)', "gradcube", norm2(x0), tf-t0, nfev(), ngev()

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if
        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("newtquad")

            call cpu_time(t0)
            call newt(x, x0, f, g, h, lsquad, gamma, eps)
            call cpu_time(tf)

            call g(x0, x)
            print '(a18, e10.2, f11.3, 2i10)', "newtquad", norm2(x0), tf-t0, nfev(), ngev()

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if
        p = fork()
        if (p == 0) then
            ! Processo filho
            call setmethod("newtcube")

            call cpu_time(t0)
            call newt(x, x0, f, g, h, lscube, gamma, eps)
            call cpu_time(tf)

            call g(x0, x)
            print '(a18, e10.2, f11.3, 3i10)', "newtcube", norm2(x0), tf-t0, nfev(), ngev(), nhev()

            deallocate(x)
            deallocate(x0)

            call exit(0)
        end if

        ! Espera os processos filhos terminarem
        do i = 1, 4
            p = wait(s)
        end do

        deallocate(x)
        deallocate(x0)
     end do
end program EP1
