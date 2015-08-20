!VFROE-ncv Approximate Riemann Solver

!SUBROUTINE VFROENCV(FR2,qR,qL,epx,epy,hmin,qs)
! Este es el que se está usando!!!!
SUBROUTINE VFROENCV(i,j,FR2,qR,qL,epR,epL,hmin,qs)
!qR=qmas
!qL=qmenos
real (kind=8), dimension(3)	:: qL,qR,qs, WL,WR,WS,WRL,Fs,WM,L,LR,LL
real (kind=8)	:: FR2, epx,epy,UL,UR,UM, hmin, alpha, beta, betaL, betaR, kappa, epxM, epyM, C
integer	::i,j, caso
real (kind=8), dimension(2)	:: epR, epL

kappa=10e-7
!L es menos
!R es mas

!Cambio de variables: 2C, u, v
WL(1)=2.0D0*sqrt(qL(1)/FR2)
WL(2)=qL(2)
WL(3)=qL(3)


WR(1)=2.0D0*sqrt(qR(1)/FR2)
WR(2)=qR(2)
WR(3)=qR(3)

! print*,'2CR=', WR(1)
! print*,'2CL=', WL(1)
! print*,'epx=', epx
! print*,'epy=', epy

!Calculo de la velocidad contravariante

! UL=WL(2)*epx+WL(3)*epy
! UR=WR(2)*epx+WR(3)*epy
! UM=0.5*(UR+UL)

UL=WL(2)*epL(1)+WL(3)*epL(2)
UR=WR(2)*epR(1)+WR(3)*epR(2)
!UM=0.5D0*(UR+UL)



!Problema Linearizado

WM(1)=0.5D0*(WL(1)+WR(1))
WM(2)=0.5D0*(WL(2)+WR(2))
WM(3)=0.5D0*(WL(3)+WR(3))

!Calculo de las metricas linearizadas
epxM=0.5D0*(epR(1)+epL(1))
epyM=0.5D0*(epR(2)+epL(2))

!Velocidad contravariante media

!UM=(WL(2)+WR(2))*0.5D0*epxM+(WL(3)+WR(3))*0.5D0*epyM

UM=WM(2)*epxM+WM(3)*epyM

!Calculo de Lambdas Linearizados
L(1)=UM-0.5D0*WM(1)*sqrt(epxM**2.0D0+epyM**2.0D0)
L(2)=UM
L(3)=UM+0.5D0*WM(1)*sqrt(epxM**2.0D0+epyM**2.0D0)

!Solver
epx=epxM
epy=epyM

WRL(:)=WR(:)-WL(:)

LR(1)=UR-0.5D0*WR(1)*sqrt(epR(1)**2.0D0+epR(2)**2.0D0)
LR(2)=UR
LR(3)=UR+0.5D0*WR(1)*sqrt(epR(1)**2.0D0+epR(2)**2.0D0)

LL(1)=UL-0.5D0*WL(1)*sqrt(epL(1)**2.0D0+epL(2)**2.0D0)
LL(2)=UL
LL(3)=UL+0.5D0*WL(1)*sqrt(epL(1)**2.0D0+epL(2)**2.0D0)

IF (qL(1)<=hmin.AND.qR(1)<=hmin) THEN

    !Caso seco
    WS(1)=0.0D0
    WS(2)=0.0D0
    WS(3)=0.0D0
    !print*, 'Caso Seco'
    caso=1

ELSE 
    !Caso clasico
  
    IF (L(1)<0.0D0.AND.L(2)<0.0D0.AND.L(3)<0.0D0) THEN	!L(1)<0, L(2)<0 y L(3)<0, W*=WR 

      WS(:)=WR(:)
      !print*, 'L3<0'
      caso=2
    ELSE IF (L(1)>0.0D0.AND.L(2)>0.0D0.AND.L(3)>0.0D0) THEN	!L(1)>0, L(2)>0 y L(3)>0, W*=WL 
      
      WS(:)=WL(:)
      caso=3
      !print*, 'L1>0'
    ELSE 

      IF (L(1)<0.0D0.AND.L(2)>0.0D0.AND.L(3)>0.0D0) THEN	!L(1)<0, L(2)>0 y L(3)>0, W*=WL+sum(L<0)L*WRL*R o sea solo L1
	
	!WS=WL+L1*WRL*R1
	
	alpha=epx**2.0D0+epy**2.0D0
	beta=-WRL(1)/(2.0D0*sqrt(alpha))+epx*WRL(2)/(2.0D0*alpha)+epy*WRL(3)/(2.0D0*alpha)
	WS(1)=WL(1)-sqrt(alpha)*beta
	WS(2)=WL(2)+epx*beta
	WS(3)=WL(3)+epy*beta
 
! 	alpha=epx**2.0D0+epy**2.0D0
! 	beta=WRL(1)/(2.0D0*sqrt(alpha))-epx*WRL(2)/(2.0D0*alpha)-epy*WRL(3)/(2.0D0*alpha)
! 	WS(1)=WL(1)+sqrt(alpha)*beta
! 	WS(2)=WL(2)-epx*beta
! 	WS(3)=WL(3)-epy*beta
	
! 	print*, 'L(1)<0, L(2)>0 y L(3)>0'	
! 	print*,'2C= ', WS(1)
	caso=41
      ELSE 		!L(1)<0, L(2)<=0 y L(3)>0, W*=WR-sum(L>0)L*WRL*R o sea solo L3
	
	!WS=WR-L3*WRL*R3
	
	alpha=epx**2.0D0+epy**2.0D0
	beta=WRL(1)/(2.0D0*sqrt(alpha))+epx*WRL(2)/(2.0D0*alpha)+epy*WRL(3)/(2.0D0*alpha)
	WS(1)=WR(1)-sqrt(alpha)*beta
	WS(2)=WR(2)-epx*beta
	WS(3)=WR(3)-epy*beta
	
 	!print*, 'ELSE'
! 	print*,'2C= ', WS(1)
	caso=42
!                IF (WS(1)<0.0D0) THEN
!		WS(1)=0
!		WS(2)=0
!		WS(3)=0
!!		print*, 'C negativo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!       print*,'2C= ', WS(1)
!	print*,'u= ', WS(2)
!	print*,'v= ', WS(3)
!	print*,'epx=',epx
!	print*,'epy=',epy
!	print*,'alpha=', alpha
!	print*,'beta=', beta
!	print*, 'WRL=', WRL
!	print*, 'WR=', WR
!	print*, 'WL=', WL
!	print*, 'FR2=', FR2
!	print*, 'qL=', qL
!	print*, 'qR=', qR
!!		write(*,'("C negativo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")') !BP
!!		read(*,*) 
!		END IF
     END IF
	
    END IF
    
    !Correcion Entropica, estado intermedio
    IF (LL(1)<0.0D0.AND.LR(1)>0.0D0) THEN
    	
! 	alpha=epx**2.0D0+epy**2.0D0
! 	betaL=-WL(1)/(2.0D0*sqrt(alpha))+WL(2)*epx/(2.0D0*alpha)+WL(3)*epy/(2.0D0*alpha)
! 	betaR=-WR(1)/(2.0D0*sqrt(alpha))+WR(2)*epx/(2.0D0*alpha)+WR(3)*epy/(2.0D0*alpha)
! 	WS(1)=WL(1)-sqrt(alpha)*0.5D0*(betaL+betaR)
! 	WS(2)=WL(2)+epx*0.5D0*(betaL+betaR)
! 	WS(3)=WL(3)+epy*0.5D0*(betaL+betaR)
 	WS(1)=WM(1)	
 	WS(2)=WM(2)	
 	WS(3)=WM(3)
	!print*,'Correccion Entropica 1'
	caso=5  
    END IF
    
    IF (LL(3)<0.0D0.AND.LR(3)>0.0D0) THEN
    	
! 	alpha=epx**2.0D0+epy**2.0D0
! 	betaL=WL(1)/(2.0D0*sqrt(alpha))+WL(2)*epx/(2.0D0*alpha)+WL(3)*epy/(2.0D0*alpha)
! 	betaR=WR(1)/(2.0D0*sqrt(alpha))+WR(2)*epx/(2.0D0*alpha)+WR(3)*epy/(2.0D0*alpha)
! 	WS(1)=WL(1)+sqrt(alpha)*0.5D0*(betaL+betaR)
! 	WS(2)=WL(2)+epx*0.5D0*(betaL+betaR)
! 	WS(3)=WL(3)+epy*0.5D0*(betaL+betaR)
	WS(1)=WM(1)	
 	WS(2)=WM(2)	
 	WS(3)=WM(3)
	!print*,'Correccion Entropica 2'
	caso=6  
    END IF

END IF

!Vuelta a variables originales

IF (WS(1)<0.0D0) THEN 
print*, 'ERROR, C is negative'
print*, i,j
	
	print*,'2C= ', WS(1)
	print*,'u= ', WS(2)
	print*,'v= ', WS(3)
	print*,'epx=',epx
	print*,'epy=',epy
	print*,'alpha=', alpha
	print*,'beta=', beta
	print*, 'WRL=', WRL
	print*, 'WR=', WR
	print*, 'WL=', WL
	print*, 'FR2=', FR2
	print*, 'qL=', qL
	print*, 'qR=', qR
	print*,'caso=',caso
!		WS(1)=0
!		WS(2)=0
!		WS(3)=0
!pause
 stop

ELSE
  IF ((0.5D0*WS(1)).le.kappa) then
  qs(1)=0.0D0 !kappa
  qs(2)=0.0D0
  qs(3)=0.0D0
  ELSE
  C=0.5D0*WS(1)
  qs(1)=((C)**2.0D0)*FR2
  qs(2)=WS(2)
  qs(3)=WS(3)
  END IF
  IF(abs(qs(2)).le.kappa) then
  qs(2)=0.0D0
  END IF
  
  IF(abs(qs(3)).le.kappa) then
  qs(3)=0.0D0
  END IF
  
END IF


END SUBROUTINE VFROENCV
