module stats
    implicit none
    save
    private
    public :: f, g, h, setf, setg, seth, initstat, printheader, printstat, setmethod, iteration, sett0, settf

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
    end interface

    procedure (f_t), pointer :: f_ptr => null()
    procedure (g_t), pointer :: g_ptr => null()
    procedure (h_t), pointer :: h_ptr => null()

    integer          :: nint
    integer          :: nf, ng, nh
    character(len=8) :: method
    double precision :: t0, tf

    double precision, allocatable :: gx(:)

contains
    subroutine initstat()
        nint = 0
        nf   = 0
        nh   = 0
        ng   = 0

        f_ptr => null()
        g_ptr => null()
        h_ptr => null()
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

    subroutine iteration()
        nint = nint + 1
    end subroutine iteration

    subroutine printheader()
        print '(18x, a16, a11, a11, a13, a14)', "‖∇f(x)‖", "t (s)", "nº f", "nº ∇f", "nº ∇²f"
    end subroutine printheader

    subroutine printstat(x)
        double precision :: x(:)

        allocate(gx(size(x)))

        call g_ptr(gx, x)
        print '(a18, e10.2, f11.3, 3i10)', method, norm2(gx), tf-t0, nf, ng, nh

        deallocate(gx)
    end subroutine printstat

    subroutine setmethod(name)
        character(len=*) :: name

        method = name
    end subroutine setmethod

    subroutine sett0()
        call cpu_time(t0)
    end subroutine sett0

    subroutine settf()
        call cpu_time(tf)
    end subroutine settf

end module stats
