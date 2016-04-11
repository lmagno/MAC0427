module entrada
contains
    subroutine le_sistema(n, A, b, filename)
        integer,           intent(out) :: n
        character(len=*),  intent(in)  :: filename
        real, allocatable, intent(out) :: b(:)
        real, allocatable, intent(out) :: A(:, :)
        
        integer :: ioerr

        ! Tenta abrir o arquivo e checa se conseguiu
        open(1, file = filename, status = 'old', iostat=ioerr)
        if (ioerr /= 0) stop "Não foi possível abrir o arquivo!"

        ! Lê o tamanho n do sistema e aloca a memória necessária
        read(1, *) n
        allocate(A(n, n))
        allocate(b(n))

        ! Lê a matriz A
        do k = 1, n*n
            read(1, *) i, j, A(i+1, j+1)
        end do

        ! Lê o vetor b
        do k = 1, n
            read(1, *) i, b(i+1)
        end do

        ! Já podemos fechar o arquivo de entrada
        close(1)
    end subroutine le_sistema
end module entrada
