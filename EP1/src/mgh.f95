module MGH
    implicit none
    save

    private
    integer :: n, nprob, ntries
    integer :: nf, ng, nh

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

    public :: setprob, getdim, gettries, getname, getinit, f, g, h, nfev, ngev, nhev
contains
    subroutine setprob(i)
        integer, intent(in)  :: i

        n      = vn(i)
        nprob  = i
        ntries = vntries(i)

        nf = 0
        ng = 0
        nh = 0
    end subroutine setprob

    function getdim()
        integer             :: getdim

        getdim = n
    end function getdim

    function gettries()
        integer :: gettries

        gettries = vntries(nprob)
    end function gettries

    function getname()
        character(len=30)   :: getname

        getname = names(nprob)
    end function getname

    subroutine getinit(x0, factor)
        double precision, intent(inout) :: x0(:)
        double precision, intent(in)    :: factor

        call INITPT(n, x0, nprob, factor)
    end subroutine getinit

    function f(x)
        double precision, intent(in) :: x(:)
        double precision             :: f

        call OBJFCN(n, x, f, nprob)
        nf = nf + 1
        if (nf > 1000000) then
            print '(i2, 3A)', nprob, "        ", names(nprob), "FC"
            call exit(1)
        end if
    end function f

    subroutine g(gx, x)
      double precision, intent(out) :: gx(:)
      double precision, intent(in)  :: x(:)

      call GRDFCN(n, x, gx, nprob)
      ng = ng + 1
  end subroutine g

    subroutine h(hx, x)
        double precision, intent(out) :: hx(:, :)
        double precision, intent(in)  :: x(:)

        integer                      :: i, j
        double precision             :: hesd(n)
        double precision             :: hesl(n*(n-1)/2)
        double precision             :: hij

        call HESFCN(n, x, hesd, hesl, nprob)

        ! Constroi a hessiana
        do j = 1, n
            hx(j, j) = hesd(j)
            do i = j+1, n
                hij     = hesl((i-1)*(i-2)/2 + j)
                hx(i, j) = hij
                hx(j, i) = hij
            end do
        end do
        nh = nh + 1
  end subroutine h

    function nfev()
        integer :: nfev

        nfev = nf
    end function nfev

    function ngev()
        integer :: ngev

        ngev = ng
    end function ngev

    function nhev()
        integer :: nhev

        nhev = nh
    end function nhev

end module MGH
