! Roda uma vez cada método (grad, newt, bfgs) com cada interpolação (lsquad, lscube)
! para o problema em questão.
! Este arquivo deve ser incluído no meio do código com a expressão INCLUDE,
! funcionando de forma similar a um macro, economizando código.

! Cria novo processo
p = fork()
if (p == 0) then
    ! Processo filho
    call setmethod("gradquad")

    call grad(x, x0, f, g, lsquad, gamma, eps, iteration, armijo)

    call printstat(x)

    deallocate(x)
    deallocate(x0)

    call exit(0)
end if

p = fork()
if (p == 0) then
    ! Processo filho
    call setmethod("gradcube")

    call grad(x, x0, f, g, lscube, gamma, eps, iteration, armijo)

    call printstat(x)

    deallocate(x)
    deallocate(x0)

    call exit(0)
end if

p = fork()
if (p == 0) then
    ! Processo filho
    call setmethod("newtquad")

    call newt(x, x0, f, g, h, lsquad, gamma, sigma, theta, eps, iteration, armijo, norm, angle)

    call printstat(x)


    deallocate(x)
    deallocate(x0)

    call exit(0)
end if

p = fork()
if (p == 0) then
    ! Processo filho
    call setmethod("newtcube")

    call newt(x, x0, f, g, h, lscube, gamma, sigma, theta, eps, iteration, armijo, norm, angle)

    call printstat(x)

    deallocate(x)
    deallocate(x0)

    call exit(0)
end if

p = fork()
if (p == 0) then
    ! Processo filho
    call setmethod("bfgsquad")

    call bfgs(x, x0, f, g, lsquad, gamma, sigma, theta, eps, iteration, armijo, norm, angle)

    call printstat(x)


    deallocate(x)
    deallocate(x0)

    call exit(0)
end if

p = fork()
if (p == 0) then
    ! Processo filho
    call setmethod("bfgscube")

    call bfgs(x, x0, f, g, lscube, gamma, sigma, theta, eps, iteration, armijo, norm, angle)

    call printstat(x)

    deallocate(x)
    deallocate(x0)

    call exit(0)
end if

! Espera os processos filhos terminarem
do i = 1, 6
    p = wait(s)
end do
