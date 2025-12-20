module dados

! Propriedades físicas do problema
real*8  :: Re  ! Número de Reynolds

! Variáveis do modelo numérico
integer :: N     ! Número de volumes de controle
integer :: iter	 ! Contador de iterações
integer :: itmax ! Número máximo de iterações
character(len=20) :: program_version = "v1.0.0-romulo" ! Versão do programa

real*8  :: dx ! Comprimento de cada volume de controle

real*8,dimension(:),allocatable :: x  ! Posição do centro de cada volume P 

real*8,dimension(:),allocatable :: ae ! Coeficiente do volume a Leste (E)
real*8,dimension(:),allocatable :: aw ! Coeficiente do volume a Oeste (W)
real*8,dimension(:),allocatable :: ap ! Coeficiente do volume central (P)
real*8,dimension(:),allocatable :: bp ! Termo-fonte do volume central (P)


! Variáveis de interesse para análise
real*8 :: um    ! Velocidade média numérica
real*8 :: uman  ! Velocidade média analítica

real*8,dimension(:),allocatable :: u    ! Velocidade no volume P (numérica)
real*8,dimension(:),allocatable :: uanl ! Velocidade no volume P (analítica)


! Variáveis do processo numérico
real*8 :: beta   ! Coeficiente do termo difusivo da função de interpolação (face w)


! Norma L1 para avaliação de erro
real*8 :: L1   ! Norma L1 do erro


! Variável para execução de comandos do sistema
integer :: dos


contains

!--------------------------------------------------------------------
! Sub-rotina para leitura dos dados de entrada
subroutine le_dados

	use portlib

	dos = system("notepad dados.txt") ! Abre o arquivo de dados no notepad (Windows)

	open(5,file='dados.txt')

	read(5,*) N
	read(5,*) Re
	read(5,*) beta
	read(5,*) itmax

	close(5)

end subroutine le_dados

!--------------------------------------------------------------------
! Sub-rotina para escrita dos dados de saída
subroutine escreve_dados

	use portlib

	integer :: i ! Contador de loop

	open(5,file='saida.txt')
	open(7,file='u.dat')

	! Cabeçalho e dados iniciais
	write(5,10) program_version
	write(5,11)
	write(5,12) N
	write(5,13) Re
	write(5,14) beta
	write(5,15) itmax

	! Coeficientes e termos-fonte calculados
	write(5,51)
	do i = 1, N
		write(5,52) i, x(i), aw(i), ae(i), ap(i), bp(i)
	end do

	! Resultados de velocidade (numérica e analítica)
	write(5,21)
	write(5,22) 1, 0.0d0, 0.0d0, 0.0d0, 0.0d0
	write(7,23) 0.0d0, 0.0d0, 0.0d0
	do i = 2, N-1
		write(5,22) i, x(i), u(i), uanl(i), uanl(i)-u(i)
		write(7,23) x(i), u(i), uanl(i)
	end do
	write(5,22) N, 1.0d0, 1.0d0, 1.0d0, 0.0d0
	write(7,23) 1.0d0, 1.0d0, 1.0d0

	! Variáveis secundárias (velocidades médias)
	write(5,34) um
	write(5,35) uman
	write(5,36) uman - um

	! Normas de erro
	write(5,41) L1
	write(5,42)	L1 / (N-2.0d0)

	close(5)
	close(7)

	! Execução de scripts de plotagem e abertura do arquivo de saída
	dos = system ("wgnuplot ucentral.gnu")
	dos = system ("wgnuplot u.gnu")
	dos = system ("notepad saida.txt")

	10 format ( "Versão do Programa: ", A )
	11 format (    "Solução Numérica do Escoamento Unidimensional em Regime Permanente:", &
	            /, "Equação de Burgers", &
	           2/, "Dados de entrada")
	12 format (/, t10, i9, t20,  " = N:      Número de volumes de controle (reais + dois fictícios)")
	13 format (  1pe18.10, t20,  " = Re:     Número de Reynolds")
	14 format (  1pe18.10, t20,  " = beta :  Coeficiente do termo difusivo (funções de interpolação):  0.0 = UDS; 1.0 = CDS")
	15 format (/, t10, i9, t20,  " = itmax:  Número máximo de iterações")

	51 format (2/, "Coeficientes e termos-fontes", &
	           2/, t3, "Volume i", t15, "xp(i)", t40, "aw(i)", t60, "ae(i)", t80, "ap(i)", t100, "bp(i)" )
	52 format (i10, 4(1pe25.15) )

	21 format (2/, "Soluções numéricas", &
	           2/, t5, "Volume", t16,"Posição [m]", t35, "Vel. numérica [m/s]", t60, "Vel. analítica [m/s]", t85, "Erro numérico [m/s]")
	22 format (i10, 1pe20.9, 3(1pe25.15) )

	23 format (1pe20.9, 2(1pe25.15) )

	34 format (2/, "Velocidade média numérica  [m/s]:      ", 1pe25.15)
	35 format (    "Velocidade média analítica [m/s]:      ", 1pe25.15)
	36 format (    "Erro numérico para a velocidade média: ", 1pe25.15)

	41 format (2/, "Norma L1:     ", 1pe25.15)
	42 format (    "Norma L1 / N: ", 1pe25.15)

end subroutine escreve_dados

!--------------------------------------------------------------------

end module dados
