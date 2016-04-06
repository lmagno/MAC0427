program EP1
    use BuscaLinear, only: maximaDescida
    use MGH,         only: setprob, getdim, gettries, getname, getinit, f, g, nfev, ngev
    implicit none

    integer                       :: p, s
    integer                       :: n, nprob, ntry, ntries
    character(len=30)             :: name
    double precision              :: factor = 1.d0
    double precision, allocatable :: x(:), x0(:)
    double precision              :: gamma, eps

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
    ! eps   = epsilon(0.d0)
    eps = 1e-4
    print '(3A)', "nprob     ", "name                          ", "‖∇f(x₀)‖  "
    do nprob = 1, 18
        ! Cria novo processo
        p = fork()
        if (p == 0) then
           call setprob(nprob)
           ntries = gettries()
           factor = 1.d0
           
           do ntry = 1, 1
              name = getname()
              n    = getdim()
              x0   = getinit(factor)
              allocate(x(n))
              
              call maximaDescida(x, x0, f, g, gamma, eps)
              print '(i2, 2A, e8.2, 2i10)', nprob, "        ", name, norm2(g(x)), nfev(), ngev()

              deallocate(x)
              deallocate(x0)
              
              factor = 10*factor
           end do

           call exit(0)
        else
           ! Espera o processo filho terminar
           p = wait(s)
        end if
     end do
end program EP1
