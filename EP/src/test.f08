module test
implicit none
contains
    function testf1(x)
        double precision, intent(in) :: x(:)
        double precision             :: testf1

        testf1 = sum(x**2)
    end function testf1

    subroutine testg1(gx, x)
        double precision, intent(out) :: gx(:)
        double precision, intent(in)  :: x(:)

        gx(:) = 2*x(:)
    end subroutine testg1

    subroutine testh1(hx, x)
        double precision, intent(out) :: hx(:, :)
        double precision, intent(in)  :: x(:)

        integer :: i, n
        n = size(x)

        hx(:, :) = 0.d0
        do i = 1, n
            hx(i, i) = 2
        end do
    end subroutine testh1

    ! x1 + x2 - 1 = 0
    subroutine c1(cx, x)
        double precision, intent(out) :: cx(:)
        double precision, intent(in)  :: x(:)

        cx(1) = x(1) + x(2) - 1
    end subroutine c1

    subroutine dc1(dcx, x)
        double precision, intent(out) :: dcx(:, :)
        double precision, intent(in)  :: x(:)

        dcx(1, 1) = 1
        dcx(1, 2) = 1
    end subroutine dc1

    subroutine d2c1(d2cx, x)
        double precision, intent(out) :: d2cx(:, :, :)
        double precision, intent(in)  :: x(:)

        d2cx(1, 1, 1) = 0
        d2cx(1, 1, 2) = 0
        d2cx(1, 2, 1) = 0
        d2cx(1, 2, 2) = 0
    end subroutine d2c1

    function testf2(x)
        double precision, intent(in) :: x(:)
        double precision             :: testf2

        testf2 = sum(x**2)
    end function testf2

    subroutine testg2(gx, x)
        double precision, intent(out) :: gx(:)
        double precision, intent(in)  :: x(:)

        gx(:) = 2*x(:)
    end subroutine testg2

    subroutine testh2(hx, x)
        double precision, intent(out) :: hx(:, :)
        double precision, intent(in)  :: x(:)

        integer :: i, n
        n = size(x)

        hx(:, :) = 0.d0
        do i = 1, n
            hx(i, i) = 2
        end do
    end subroutine testh2

    ! x1 + x2 - 1 = 0
    subroutine c2(cx, x)
        double precision, intent(out) :: cx(:)
        double precision, intent(in)  :: x(:)

        cx(1) = (x(1) + x(2) - 1)**3
    end subroutine c2

    subroutine dc2(dcx, x)
        double precision, intent(out) :: dcx(:, :)
        double precision, intent(in)  :: x(:)

        dcx(1, 1) = 3*(x(1) + x(2) - 1)**2
        dcx(1, 2) = 3*(x(1) + x(2) - 1)**2
    end subroutine dc2

    subroutine d2c2(d2cx, x)
        double precision, intent(out) :: d2cx(:, :, :)
        double precision, intent(in)  :: x(:)

        d2cx(1, 1, 1) = 6*(x(1) + x(2) - 1)
        d2cx(1, 1, 2) = 6*(x(1) + x(2) - 1)
        d2cx(1, 2, 1) = 6*(x(1) + x(2) - 1)
        d2cx(1, 2, 2) = 6*(x(1) + x(2) - 1)
    end subroutine d2c2

    function testf3(x)
        double precision, intent(in) :: x(:)
        double precision             :: testf3

        testf3 = sum(x**2)
    end function testf3

    subroutine testg3(gx, x)
        double precision, intent(out) :: gx(:)
        double precision, intent(in)  :: x(:)

        gx(:) = 2*x(:)
    end subroutine testg3

    subroutine testh3(hx, x)
        double precision, intent(out) :: hx(:, :)
        double precision, intent(in)  :: x(:)

        integer :: i, n
        n = size(x)

        hx(:, :) = 0.d0
        do i = 1, n
            hx(i, i) = 2
        end do
    end subroutine testh3

    ! x1 + x2 - 1 = 0
    subroutine c3(cx, x)
        double precision, intent(out) :: cx(:)
        double precision, intent(in)  :: x(:)

        cx(1) = x(1) - 1
        cx(2) = x(2)
    end subroutine c3

    subroutine dc3(dcx, x)
        double precision, intent(out) :: dcx(:, :)
        double precision, intent(in)  :: x(:)

        dcx(1, 1) = 1
        dcx(1, 2) = 0

        dcx(2, 1) = 0
        dcx(2, 2) = 1
    end subroutine dc3

    subroutine d2c3(d2cx, x)
        double precision, intent(out) :: d2cx(:, :, :)
        double precision, intent(in)  :: x(:)

        d2cx(1, 1, 1) = 0
        d2cx(1, 1, 2) = 0
        d2cx(1, 2, 1) = 0
        d2cx(1, 2, 2) = 0

        d2cx(2, 1, 1) = 0
        d2cx(2, 1, 2) = 0
        d2cx(2, 2, 1) = 0
        d2cx(2, 2, 2) = 0
    end subroutine d2c3

    function testf4(x)
        double precision, intent(in) :: x(:)
        double precision             :: testf4

        testf4 = sum(x**2)
    end function testf4

    subroutine testg4(gx, x)
        double precision, intent(out) :: gx(:)
        double precision, intent(in)  :: x(:)

        gx(:) = 2*x(:)
    end subroutine testg4

    subroutine testh4(hx, x)
        double precision, intent(out) :: hx(:, :)
        double precision, intent(in)  :: x(:)

        integer :: i, n
        n = size(x)

        hx(:, :) = 0.d0
        do i = 1, n
            hx(i, i) = 2
        end do
    end subroutine testh4

    ! x1 + x2 - 1 = 0
    subroutine c4(cx, x)
        double precision, intent(out) :: cx(:)
        double precision, intent(in)  :: x(:)

        cx(1) = (x(1) - 1)**2 + (x(2) - 1)**2 - 1
        cx(2) = (x(1) - 1)**2 + (x(2) + 1)**2 - 1
    end subroutine c4

    subroutine dc4(dcx, x)
        double precision, intent(out) :: dcx(:, :)
        double precision, intent(in)  :: x(:)

        dcx(1, 1) = 2*(x(1) - 1)
        dcx(1, 2) = 2*(x(2) - 1)

        dcx(2, 1) = 2*(x(1) - 1)
        dcx(2, 2) = 2*(x(2) + 1)
    end subroutine dc4

    subroutine d2c4(d2cx, x)
        double precision, intent(out) :: d2cx(:, :, :)
        double precision, intent(in)  :: x(:)

        d2cx(1, 1, 1) = 2
        d2cx(1, 1, 2) = 0
        d2cx(1, 2, 1) = 0
        d2cx(1, 2, 2) = 2

        d2cx(2, 1, 1) = 2
        d2cx(2, 1, 2) = 0
        d2cx(2, 2, 1) = 0
        d2cx(2, 2, 2) = 2
    end subroutine d2c4


end module test
