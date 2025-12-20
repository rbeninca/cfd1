module coeficientes

use dados

contains

!--------------------------------------------------------------------
! Sub-rotina para definir coeficientes fictícios de contorno
subroutine coef_fict

	! Coeficientes para o primeiro volume (condição de contorno de entrada)
	aw(1) =  0.0d0
	ae(1) = -1.0d0
	ap(1) =  1.0d0
	bp(1) =  0.0d0

	! Coeficientes para o último volume (condição de contorno de saída)
	aw(N) = -1.0d0
	ae(N) =  0.0d0
	ap(N) =  1.0d0
	bp(N) =  2.0d0

end subroutine coef_fict

!--------------------------------------------------------------------
! Sub-rotina para calcular coeficientes reais para volumes internos
subroutine coef_reais

	integer :: i   ! Variável de iteração
	real*8  :: S   ! Variável auxiliar para termo fonte
    real*8  :: local_dx ! Cópia local de dx para demonstração

    local_dx = dx ! Atribuição de dx a uma variável local

	! Loop sobre os volumes internos
	do i = 2, N-1
		aw(i) = 2.0d0 + Re*( u(i-1) + u(i) )*local_dx
		ae(i) = 2.0d0
		ap(i) =	4.0d0 + Re*( u(i) + u(i+1) )*local_dx

		S = (Re**2)*dexp( x(i)*Re )*( 2.0d0*dexp( x(i)*Re )-dexp(Re)-1.0d0 ) / &
		    ( dexp(Re)-1.0d0 )**2 
		bp(i) = 2.0d0*S*(local_dx**2) + 0.5d0*beta*Re*( 2.0d0*u(i)**2 - u(i-1)**2 - u(i+1)**2 )*local_dx
	end do

end subroutine coef_reais

!--------------------------------------------------------------------
! Sub-rotina para definir as posições dos volumes de controle
subroutine posicao

	integer :: i ! Variável de iteração

	x(1) = 0.0d0
	do i = 2, N-1
		x(i) = (i-1.5d0)*dx
	end do
	x(N) = 1.0d0

end subroutine posicao

!--------------------------------------------------------------------

end module coeficientes
