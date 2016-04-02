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

    public :: setprob, getdim, getname, getinit, f, g
contains
    subroutine setprob(i)
        integer, intent(in)  :: i

        n      = vn(i)
        nprob  = i
        ntries = vntries(i)
    end subroutine setprob

    function getdim()
        integer             :: getdim

        getdim = n
    end function getdim

    function getname()
        character(len=30)   :: getname

        getname = names(nprob)
    end function getname

    function getinit(factor)
        double precision, intent(in)  :: factor
        double precision              :: getinit(n)

        call INITPT(n, getinit, nprob, factor)
    end function getinit

    function f(x)
        double precision, intent(in) :: x(:)
        double precision             :: f

        call OBJFCN(n, x, f, nprob)
    end function f

    function g(x)
      double precision, intent(in) :: x(:)
      double precision             :: g(size(x))

      call GRDFCN(n, x, g, nprob)
    end function g
end module MGH
