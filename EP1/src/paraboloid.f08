module paraboloid
    implicit none

contains
    function paraf(x)
        double precision, intent(in) :: x(:)
        double precision             :: paraf

        double precision :: x1, x2
        double precision :: a, b

        x1 = x(1)
        x2 = x(2)

        paraf = x1*x1 + x2*x2
    end function paraf

    subroutine parag(gx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: gx(:)

        double precision :: x1, x2

        x1 = x(1)
        x2 = x(2)

        gx(1) = 2*x1
        gx(2) = 2*x2
    end subroutine parag

    subroutine parah(hx, x)
        double precision, intent(in)  :: x(:)
        double precision, intent(out) :: hx(:, :)

        double precision :: x1, x2

        x1 = x(1)
        x2 = x(2)

        hx(1, 1) = 2
        hx(1, 2) = 0
        hx(2, 1) = 0
        hx(2, 2) = 2
    end subroutine parah
end module paraboloid
