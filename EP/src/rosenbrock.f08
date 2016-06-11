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
