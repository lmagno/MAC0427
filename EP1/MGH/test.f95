program main
    implicit none

    integer :: n, m, lfj, mode
    real    :: ftf
    real, allocatable    :: x(:), f(:), fj(:,:), g(:)

    integer            nprob, nprobs, nstart, nstrts
    common /PROBLM/    nprob, nprobs, nstart, nstrts

    integer            nout
    common /IOUNIT/    nout

    real               coeff
    common /PARAM1/    coeff

    mode = -1
    nout = 6

    call getfun(x, n, f, m, ftf, fj, lfj, g, mode)

    print *, nprob, nprobs, nstart, nstrts
    print *, nout
    print *, coeff
end program main
