module constrained
    use       ls, only: lsquad
    use   trisys, only: backcol, forwcol
    use cholesky, only: cholcol
implicit none
    integer :: n           ! Dimensão do problema
    integer :: m           ! Número de restrições do problema
    double precision :: mu ! Peso da penalidade

    interface
        function f_t(x)
            double precision, intent(in) :: x(:)
            double precision             :: f_t
        end function f_t

        subroutine c_t(cx, x)
            double precision, intent(out) :: cx(:)
            double precision, intent(in)  :: x(:)
        end subroutine c_t
    end interface

    procedure (f_t), pointer :: f_ptr => null()
    procedure (c_t), pointer :: c_ptr => null()
contains

    function Q(x)
        double precision, intent(in) :: x(:)
        double precision             :: Q

        double precision :: cx(size(x))

        call c_ptr(cx, x)
        Q = f_ptr(x) + 1/(2*mu)*sum(cx**2)
    end function Q

    subroutine dQ(dQx, gx, cx, dcx)
        double precision, intent(out) :: dQx(:)
        double precision, intent(in)  :: gx(:)
        double precision, intent(in)  :: cx(:)
        double precision, intent(in)  :: dcx(:, :)
        integer :: k, i

        dQx(:) = gx(:)

        do k = 1, m
            do i = 1, n
                dQx(i) = dQx(i) + cx(k)*dcx(k, i)/mu
            end do
        end do
    end subroutine dQ

    subroutine S(Sx, hx, cx, dcx, d2cx)
        double precision, intent(out) :: Sx(:, :)
        double precision, intent(in)  :: hx(:, :)
        double precision, intent(in)  :: cx(:)
        double precision, intent(in)  :: dcx(:, :)
        double precision, intent(in)  :: d2cx(:, :, :)
        integer :: i, j, k

        do j = 1, n
            do i = 1, n
                Sx(i, j) = hx(i, j)

                do k = 1, m
                    Sx(i, j) = Sx(i, j) + cx(k)*d2cx(k, i, j)/mu
                end do
            end do

            do i = n+1, n+m
                Sx(i, j) = dcx(i-n, j)
            end do
        end do

        do j = n+1, n+m
            do i = 1, n
                Sx(i, j) = dcx(j-n, i)
            end do

            do i = n+1, n+m
                if (i == j) then
                    Sx(i, j) = -mu
                else
                    Sx(i, j) = 0.d0
                end if
            end do
        end do
    end subroutine S

    subroutine v(vx, dQx)
        double precision, intent(out) :: vx(:)
        double precision, intent(in) :: dQx(:)
        integer i

        do i = 1, n
            vx(i) = -dQx(i)
        end do

        do i = n+1, n+m
            vx(i) = 0.d0
        end do
    end subroutine v

    subroutine newt(x, f, g, h, c, dc, d2c, gamma, sigma, theta, eps, iteration, armijo, norm, angle)
        ! Saída
        double precision, intent(inout) :: x(:)

        ! Entrada
        double precision, intent(in) :: eps     ! Tolerância
        double precision, intent(in) :: gamma   ! γ > 0, Constante para condição de Armijo
        double precision, intent(in) :: sigma   ! σ > 0, constante para condição da norma
        double precision, intent(in) :: theta   ! Θ > 0, constante para condição do ângulo

        ! Função, gradiente, hessiana e busca linear
        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine g(gx, x)
                double precision, intent(out) :: gx(:)
                double precision, intent(in)  :: x(:)
            end subroutine g

            subroutine h(hx, x)
                double precision, intent(out) :: hx(:, :)
                double precision, intent(in)  :: x(:)
            end subroutine h

            subroutine c(cx, x)
                double precision, intent(out) :: cx(:)
                double precision, intent(in)  :: x(:)
            end subroutine c

            subroutine dc(dcx, x)
                double precision, intent(out) :: dcx(:, :)
                double precision, intent(in)  :: x(:)
            end subroutine dc

            subroutine d2c(d2cx, x)
                double precision, intent(out) :: d2cx(:, :, :)
                double precision, intent(in)  :: x(:)
            end subroutine d2c

            subroutine iteration()
            end subroutine iteration

            subroutine armijo()
            end subroutine armijo

            subroutine norm()
            end subroutine norm

            subroutine angle()
            end subroutine angle
        end interface


        ! Variáveis locais
        double precision :: d(n)    ! Direção
        double precision :: dz(n+m) ! Direção aumentada
        integer          :: p(n)    ! Vetor de permutação
        double precision :: alpha   ! Passo
        double precision :: rho     ! Constante para shift nos autovalores
                                    ! da hessiana

        integer          :: i, status, tmp

        ! Variáveis para evitar chamadas de função desnecessárias
        double precision :: xd(n)
        double precision :: Sr(n+m, n+m) ! S(x)
        double precision :: dQTd         ! ∇ᵀf(x)d
        double precision :: ndQ          ! ‖∇f(x)‖
        double precision :: nd           ! ‖d‖

        double precision :: fx            ! f(x)
        double precision :: gx(n)         ! ∇f(x)
        double precision :: hx(n, n)      ! ∇²f(x)

        double precision :: cx(m)         ! [c₁(x) ..   cₘ(x)]ᵀ
        double precision :: dcx(m, n)     ! [∇c₁(x) ..  ∇cₘ(x)]ᵀ = A(x)
        double precision :: d2cx(m, n, n) ! [∇²c₁(x) .. ∇²cₘ(x)]ᵀ

        double precision :: Qx            !   Q(x, μ) =     f(x) + 1/2μ*∑cᵢ(x)²
        double precision :: dQx(n)        !  ∇Q(x, μ) =   ∇f(x) + 1/μ*∑cᵢ(x)∇cᵢ(x)
        double precision :: d2Qx(n, n)    ! ∇²Q(x, μ) =  ∇²f(x) + 1/μ*∑cᵢ(x)∇²cᵢ(x) + 1/μ*A(x)ᵀA(x)

        double precision :: vx(n+m), Sx(n+m, n+m) ! S(x)*[d z]ᵀ = v(x)


        ! Valores iniciais
        call g(gx, x)
        call c(cx, x)
        call dc(dcx, x)
        call dQ(dQx, gx, cx, dcx)
        call v(vx, dQx)

        do while (norm2(dQx) > eps)
            call iteration()

            ! Calcula fator de Cholesky de uma aproximação
            ! definida positiva de S(x)
            !     GGᵀ = S(x) + ρI
            call h(hx, x)
            call d2c(d2cx, x)
            call S(Sx, hx, cx, dcx, d2cx)

            ! Verifica se a própria S(x) não é definida positiva
            Sr  = Sx
            status = cholcol(n+m, Sr)
            if (status == 0) then
                ! É definida positiva

                dz = vx
                ! Resolve Gᵀdz solução da equação
                !    GGᵀdz = v(x)
                tmp = forwcol(n+m, Sr, dz, unit = .false.)

                ! Resolve dz solução da equação
                !     Gᵀdz = G⁻¹v(x)
                tmp = backcol(n+m, Sr, dz, trans = .true.)
            end if

            ! Extrai a direção de descida do vetor aumentado
            do i = 1, n
                d(i) = dz(i)
            end do

            ndQ  = norm2(dQx)
            nd  = norm2(d)
            dQTd = dot_product(dQx, d)
            rho = 1e-3

            ! Condição do ângulo
            do while (status == -1 .or. dQTd > -theta*ndQ*nd)
                call angle()
                ! Faz shift nos autovalores e tenta novamente
                ! Sr = S(x) + ρI
                Sr = Sx
                do i = 1, n+m
                    Sr(i, i) = Sr(i, i) + rho
                end do

                status = cholcol(n+m, Sr)
                if (status == 0) then
                    ! É definida positiva

                    dz = vx
                    ! Resolve Gᵀd solução da equação
                    !    GGᵀdz = v(x)
                    tmp = forwcol(n+m, Sr, dz, unit = .false.)

                    ! Resolve d solução da equação
                    !     Gᵀdz = G⁻¹v(x)
                    tmp = backcol(n+m, Sr, dz, trans = .true.)
                end if

                ! Extrai a direção de descida do vetor aumentado
                do i = 1, n
                    d(i) = dz(i)
                end do

                ndQ  = norm2(dQx)
                nd   = norm2(d)
                dQTd = dot_product(dQx, d)
                rho  = 10*rho
            end do

            ! Condição da norma
            if (nd < sigma*ndQ) then
                call norm()
                d = sigma*ndQ/nd*d
            end if

            dQTd = dot_product(dQx, d)

            ! Busca linear
            alpha = 1.d0
            call lsquad(alpha, x, d, Q, gamma, dQTd, xd, armijo)

            ! Atualiza x com o novo valor
            x  = x + alpha*d
            call g(gx, x)
            call c(cx, x)
            call dc(dcx, x)
            call dQ(dQx, gx, cx, dcx)
            call v(vx, dQx)
        end do
    end subroutine newt

    subroutine penalty(x, x0, n_, m_, mu_, f, g, h, c, dc, d2c, &
                        gamma, sigma, theta, eps, subproblem, iteration, armijo, norm, angle)
        ! Entrada
        double precision, intent(in)    :: x0(:)   ! Ponto inicial
        integer,          intent(in)    :: n_, m_
        double precision, intent(inout) :: mu_     ! μ > 0, valor inicial para o peso da penalidade
        double precision, intent(in)    :: gamma   ! γ > 0, Constante para condição de Armijo
        double precision, intent(in)    :: sigma   ! σ > 0, constante para condição da norma
        double precision, intent(in)    :: theta   ! Θ > 0, constante para condição do ângulo
        double precision, intent(in)    :: eps     ! Tolerância

        ! Função, gradiente, hessiana e busca linear
        interface
            function f(x)
                double precision, intent(in) :: x(:)
                double precision             :: f
            end function f

            subroutine g(gx, x)
                double precision, intent(out) :: gx(:)
                double precision, intent(in)  :: x(:)
            end subroutine g

            subroutine h(hx, x)
                double precision, intent(out) :: hx(:, :)
                double precision, intent(in)  :: x(:)
            end subroutine h

            subroutine c(cx, x)
                double precision, intent(out) :: cx(:)
                double precision, intent(in)  :: x(:)
            end subroutine c

            subroutine dc(dcx, x)
                double precision, intent(out) :: dcx(:, :)
                double precision, intent(in)  :: x(:)
            end subroutine dc

            subroutine d2c(d2cx, x)
                double precision, intent(out) :: d2cx(:, :, :)
                double precision, intent(in)  :: x(:)
            end subroutine d2c

            subroutine subproblem()
            end subroutine subproblem

            subroutine iteration()
            end subroutine iteration

            subroutine armijo()
            end subroutine armijo

            subroutine norm()
            end subroutine norm

            subroutine angle()
            end subroutine angle
        end interface

        ! Saída
        double precision, intent(out) :: x(:)

        ! Variáveis locais
        double precision :: cx(m_)
        double precision :: dcx(m_, n_)
        double precision :: gx(n_)
        double precision :: dL(n_) ! ∇L(x, -c(x)/μ) = ∇f(x) + 1/μ*∑cₖ(x)∇cₖ(x)
        integer          :: i, j, k
        logical          :: kkt

        n = n_
        m = m_
        mu = mu_
        f_ptr => f
        c_ptr => c

        ! Valores iniciais
        x(:) = x0(:)
        call c(cx, x)
        call dc(dcx, x)
        call g(gx, x)

        ! Lagrangeana
        dL(:) = gx(:)
        do k = 1, m
            do i = 1, n
                dL(i) = dL(i) + cx(k)*dcx(k, i)/mu
            end do
        end do

        write (19, *), x
        if (norm2(cx) < eps .and. norm2(dL) < eps) then
            kkt = .true.
        else
            kkt = .false.
        end if

        do while(.not. kkt)
            ! Minimiza Q(x, μ)
            call subproblem()
            call newt(x, f, g, h, c, dc, d2c, gamma, sigma, theta, eps, iteration, armijo, norm, angle)

            ! Verifica KKT
            call c(cx, x)
            call dc(dcx, x)
            call g(gx, x)

            dL(:) = gx(:)
            do k = 1, m
                do i = 1, n
                    dL(i) = dL(i) + cx(k)*dcx(k, i)/mu
                end do
            end do

            write (19, *), x
            write (38, '(e10.2, a, e10.2, a)'), norm2(cx), " & ", norm2(dL), " \\"
            if (norm2(cx) < eps .and. norm2(dL) < eps) then
                kkt = .true.
            else
                ! Atualiza μ
                mu = mu/2
            end if
        end do

        mu_ = mu
    end subroutine penalty

end module constrained
