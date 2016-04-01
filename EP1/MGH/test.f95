module MGH
    implicit none
    save

    private
    integer :: n, nprob, ntries

    integer :: nprobs  = 18
    integer,           dimension(18) :: vn      = (/ 3, 6, 3, 2, 3, 4, 10, 2, 2, 2, 4, 3, 2, 12, 4, 2, 4, 2/)
    integer,           dimension(18) :: vntries = (/ 2, 2, 2, 3, 2, 4,  4, 2, 2, 5, 2, 2, 2,  5, 2, 2, 4, 2/)
    character(len=30), dimension(18) :: names   = [character(len=30) :: &
                                                   "Helical Valley", &
                                                   "Biggs EXP6", &
                                                   "Gaussian", &
                                                   "Powell Badly Scaled", &
                                                   "Box 3D", &
                                                   "Variably Dimensioned", &
                                                   "Watson", &
                                                   "Penalty I", &
                                                   "Penalty II", &
                                                   "Brown Badly Scaled", &
                                                   "Brown and Dennis", &
                                                   "Gulf Research and Development", &
                                                   "Trigonometric", &
                                                   "Extended Rosebrock", &
                                                   "Extended Powell Singular", &
                                                   "Beale", &
                                                   "Wood", &
                                                   "Chebyquad"]

    public :: setprob, getdim, getname, setinit, f
contains
    subroutine setprob(mprob)
        integer, intent(in)  :: mprob

        n      = vn(mprob)
        nprob  = mprob
        ntries = vntries(mprob)
    end subroutine setprob

    function getdim()
        integer             :: getdim

        getdim = n
    end function getdim

    function getname()
        character(len=30)   :: getname

        getname = names(nprob)
    end function getname

    subroutine setinit(x, factor)
        double precision, intent(in)  :: factor
        double precision, intent(out) :: x(n)

        call INITPT(n, x, nprob, factor)
    end subroutine setinit

    function f(x)
        double precision :: x(n)
        double precision :: f

        call OBJFCN(n, x, f, nprob)
    end function f
end module MGH

program main
    use MGH
    implicit none

    integer                       :: n, nprob
    character(len=30)             :: name
    double precision              :: factor = 1.d0
    double precision, allocatable :: x(:)

    call setprob(3)
    name = getname()
    print *, name
    n = getdim()
    print *, "n =", n
    allocate(x(n))
    x(:) = 0.d0
    print *, "x     =", x
    print *, "f(x)  =", f(x)

    call setinit(x, factor)
    print *, "x0    =", x
    print *, "f(x0) =", f(x)

    ! print *, vn
    ! print *, names(1)
    deallocate(x)
end program main
