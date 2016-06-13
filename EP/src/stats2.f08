module stats2
    implicit none
    save
    private
    public :: f, g, h, &
              setf, setg, seth, &
              setc, setdc, setd2c, &
              iteration, armijo, norm, angle, &
              initstat, printheader, printstat, setmethod

    interface
        function f_t(x)
            double precision, intent(in) :: x(:)
            double precision             :: f_t
        end function f_t

        subroutine g_t(gx, x)
            double precision, intent(in)  :: x(:)
            double precision, intent(out) :: gx(:)
        end subroutine g_t

        subroutine h_t(hx, x)
            double precision, intent(in)  :: x(:)
            double precision, intent(out) :: hx(:, :)
        end subroutine h_t


        subroutine c_t(cx, x)
            double precision, intent(in)  :: x(:)
            double precision, intent(out) :: cx(:)
        end subroutine c_t

        subroutine dc_t(dcx, x)
            double precision, intent(in)  :: x(:)
            double precision, intent(out) :: dcx(:, :)
        end subroutine dc_t

        subroutine d2c_t(d2cx, x)
            double precision, intent(in)  :: x(:)
            double precision, intent(out) :: d2cx(:, :, :)
        end subroutine d2c_t
    end interface

    procedure   (f_t), pointer ::   f_ptr => null()
    procedure   (g_t), pointer ::   g_ptr => null()
    procedure   (h_t), pointer ::   h_ptr => null()
    procedure   (c_t), pointer ::   c_ptr => null()
    procedure  (dc_t), pointer ::  dc_ptr => null()
    procedure (d2c_t), pointer :: d2c_ptr => null()

    integer          :: n, m
    integer          :: nit, narm, nnor, nang
    integer          :: nf, ng, nh
    integer          :: nc, ndc, nd2c
    character(len=8) :: method

    double precision, allocatable :: gx(:)

contains
    subroutine initstat(n_, m_)
        integer, intent(in) :: n_, m_

        n = n_
        m = m_
        nit = 0
        narm = 0
        nnor = 0
        nang = 0
        nf   = 0
        nh   = 0
        ng   = 0

        method = ""
          f_ptr => null()
          g_ptr => null()
          h_ptr => null()
          c_ptr => null()
         dc_ptr => null()
        d2c_ptr => null()
    end subroutine initstat

    subroutine setf(p)
        procedure (f_t) :: p
        f_ptr => p
    end subroutine setf

    subroutine setg(p)
        procedure (g_t) :: p
        g_ptr => p
    end subroutine setg

    subroutine seth(p)
        procedure (h_t) :: p
        h_ptr => p
    end subroutine seth

    subroutine setc(p)
        procedure (c_t) :: p
        c_ptr => p
    end subroutine setc

    subroutine setdc(p)
        procedure (dc_t) :: p
        dc_ptr => p
    end subroutine setdc

    subroutine setd2c(p)
        procedure (d2c_t) :: p
        d2c_ptr => p
    end subroutine setd2c

    function f(x)
        double precision, intent(in) :: x(:)
        double precision             :: f

        nf = nf + 1

        if (nf > 1000000) then
            print '(a18, a10)', method, "FC"
            call exit(1)
        end if

        f = f_ptr(x)
    end function f

    subroutine g(gx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: gx(:)

        ng = ng + 1

        call g_ptr(gx, x)
    end subroutine g

    subroutine h(hx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: hx(:, :)

        nh = nh + 1

        call h_ptr(hx, x)
    end subroutine h

    subroutine c(cx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: cx(:)

        nc = nc + 1

        call c_ptr(cx, x)
    end subroutine c

    subroutine dc(dcx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: dcx(:, :)

        ndc = ndc + 1

        call dc_ptr(dcx, x)
    end subroutine dc

    subroutine d2c(d2cx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: d2cx(:, :, :)

        nd2c = nd2c + 1

        call d2c_ptr(d2cx, x)
    end subroutine d2c

    subroutine iteration()
        nit = nit + 1
    end subroutine iteration

    subroutine armijo()
        narm = narm + 1
    end subroutine armijo

    subroutine norm()
        nnor = nnor + 1
    end subroutine norm

    subroutine angle()
        nang = nang + 1
    end subroutine angle

    subroutine printheader()
        print '(a22, a16, a10, a12, a13, 2a10, a11)', "‖c(x)‖", "‖∇L‖", "it", "f", "∇f", "∇²f", &
                                                           "armijo", "norma", "ângulo"
    end subroutine printheader

    subroutine printstat(x, mu)
        double precision, intent(in) :: x(:)
        double precision, intent(in) :: mu

        double precision, allocatable :: cx(:)
        double precision, allocatable :: dcx(:, :)
        double precision, allocatable :: gx(:)
        double precision, allocatable :: dL(:)

        integer :: i, k

        allocate(cx(m))
        allocate(dcx(m, n))
        allocate(gx(n))
        allocate(dL(n))

        call c_ptr(cx, x)
        call dc_ptr(dcx, x)
        call g_ptr(gx, x)

        dL(:) = gx(:)
        do k = 1, m
            do i = 1, n
                dL(i) = dL(i) + cx(k)*dcx(k, i)/mu
            end do
        end do
        print '(a8, 2e10.2, 7i10)', method, norm2(cx), norm2(dL), nit, nf, ng, nh, narm, nnor, nang

        deallocate(cx)
        deallocate(dcx)
        deallocate(gx)
        deallocate(dL)
    end subroutine printstat

    subroutine setmethod(name)
        character(len=*) :: name

        method = name
    end subroutine setmethod

end module stats2
