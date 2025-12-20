module resultados

use dados
use coeficientes
use solvers

contains

!--------------------------------------------------------------------
! Sub-rotina principal para o processamento da simulação
subroutine processamento

	real*8 :: ucent ! Variável auxiliar para u(x=1/2) da iteração anterior

	! Leitura dos dados de entrada
	call le_dados

	! Alocação dinâmica das variáveis
	allocate ( aw(N), ae(N), ap(N), bp(N), u(N), uanl(N), x(N) )

	! Definição do comprimento do volume de controle
	dx = 1.0d0/(N-2.0d0)

	! Definição da posição de cada volume de controle
	call posicao

	! Cálculo da solução analítica
	call solucao_analitica

	! Inicialização da solução numérica com a solução analítica
	u = uanl
	ucent = u((N+1)/2)

	! Cálculo dos coeficientes e termos-fonte fictícios (condições de contorno)
	call coef_fict

	open(8,file = 'ucentral.dat')

	! Loop principal de iterações
	do iter = 1, itmax

		! Cálculo dos coeficientes e termos-fonte reais para os volumes internos
		call coef_reais

		! Solução do sistema de equações lineares usando o método TDMA
		call tdma ( N, ap, aw, ae, bp, u )

		! Escrita do erro de convergência para o ponto central
		write(8,5) iter, dabs( ucent-u((N+1)/2) )
		ucent = u((N+1)/2)

        ! Exemplo de linha de depuração (não altera a lógica principal)
        print *, "Iteração: ", iter, " | u_central: ", u((N+1)/2)

	end do

	close (8)

	! Pós-processamento dos resultados
	call pos_processamento

	! Escrita dos dados de saída para arquivos
	call escreve_dados

	5 format ( i9, 1pe25.15 )

end subroutine processamento

!--------------------------------------------------------------------
! Sub-rotina para pós-processamento dos resultados
subroutine pos_processamento

	integer :: i ! Contador de loop

	! Cálculo da velocidade média numérica
	um = 0.0d0
	do i = 2, N-1
		um = um + u(i)
	end do
	um = um*dx 

	! Cálculo da Norma L1 do erro
	L1 = 0.0d0
	do i = 2, N-1
		L1 = L1 + dabs( uanl(i) - u(i) )
	end do

	! Cálculo da velocidade média analítica
	uman = (dexp(Re) - Re - 1.0d0) / ( Re*(dexp(Re)-1.0d0) )

end subroutine pos_processamento

!--------------------------------------------------------------------
! Sub-rotina para cálculo da solução analítica
subroutine solucao_analitica

	integer :: i ! Contador de loop
	
	! Cálculo da velocidade analítica para cada ponto
	do i = 1, N
		uanl(i) =  (dexp(x(i)*Re) - 1.0d0) / (dexp(Re) - 1.0d0)
	end do

end subroutine solucao_analitica

!--------------------------------------------------------------------

end module resultados
