SUBROUTINE friccion(q,Nx,Ny,dt,Fr2,hmin,kappa,qf,Cf,C)
! q(h,u,v)
! Recibe la matriz de coeficientes de friccion C asociado a cada nodo (i,j)
! Friction: Manning and Chezy Approach

integer :: Nx, Ny, Cf
real (kind=8),dimension(3,Nx,Ny):: q, qf, qfaux
real (kind=8):: dt, Fr2, F1, F2, alpha1, alpha2, beta1, beta2, qaux1, qaux2
real (kind=8),dimension(Nx,Ny):: C
real (kind=8):: h0tmp
!allocate(q(3,Nx,Ny),qf(3,Nx,Ny),qfaux(3,Nx,Ny))

h0tmp=0.5D0

DO i=1,Nx; DO j=1,Ny

!qaux(h,hu,hv)
qaux1=q(1,i,j)*q(2,i,j)
qaux2=q(1,i,j)*q(3,i,j)

!qaux1=q(2,i,j)
!qaux2=q(3,i,j)


IF (q(1,i,j)==0.0D0) then !Caso Seco

F1=0.0D0
F2=0.0D0

ELSE IF (q(1,i,j)/=0.0D0.AND.q(2,i,j)==0.0D0.AND.q(3,i,j)==0.0D0) then !Reposo
F1=0.0D0
F2=0.0D0

ELSE

alpha1=(1.0D0/Fr2)*q(2,i,j)*sqrt(q(2,i,j)**2.0D0+q(3,i,j)**2.0D0)
alpha2=(1.0D0/Fr2)*q(3,i,j)*sqrt(q(2,i,j)**2.0D0+q(3,i,j)**2.0D0)
beta1=(2.0D0*qaux1**2.0D0+qaux2**2.0D0)/sqrt(qaux1**2.0D0+qaux2**2.0D0)*(1/Fr2*q(1,i,j)**2.0D0)
beta2=(qaux1**2.0D0+2.0D0*qaux2**2.0D0)/sqrt(qaux1**2.0D0+qaux2**2.0D0)*(1/Fr2*q(1,i,j)**2.0D0)

    IF (Cf==1) THEN  !Manning
    !print *, h0tmp
    F1=(-(C(i,j)**2.0D0/q(1,i,j)**(1.0D0/3.0D0))*alpha1)/(1.0D0+dt*(C(i,j)**2.0D0/q(1,i,j)**(1.0D0/3.0D0))*beta1)
    F2=(-(C(i,j)**2.0D0/q(1,i,j)**(1.0D0/3.0D0))*alpha2)/(1.0D0+dt*(C(i,j)**2.0D0/q(1,i,j)**(1.0D0/3.0D0))*beta2)
!!!!!!!!!!!!!!!!!!!!!!!!! ESSAI AVEC MANNING AVES HAUTEUR CONSTANTE A CHANGER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !F1=(-(C(i,j)**2.0D0/q(1,i,j)**(1.0D0/3.0D0))*alpha1)/(1.0D0+2.0D0*dt/Fr2*(C(i,j)**2.0D0/q(1,i,j)**(4.0D0/3.0D0))*q(2,i,j))
    !F2=(-(C(i,j)**2.0D0/q(1,i,j)**(1.0D0/3.0D0))*alpha2)/(1.0D0+2.0D0*dt/Fr2*(C(i,j)**2.0D0/q(1,i,j)**(4.0D0/3.0D0))*q(3,i,j))
    !F1= -alpha1*C(i,j)*Fr2/(1.0D0+dt*C(i,j)*Fr2*beta1)
    !F2= -alpha2*C(i,j)*Fr2/(1.0D0+dt*C(i,j)*Fr2*beta2)
    ELSE IF (Cf==2) THEN !Chezy
    
    F1=(-(1.0D0/C(i,j)**2.0D0)*alpha1)/(1.0D0+dt*(1.0D0/C(i,j)**2.0D0)*beta1)
    F2=(-(1.0D0/C(i,j)**2.0D0)*alpha2)/(1.0D0+dt*(1.0D0/C(i,j)**2.0D0)*beta2)

    ELSE IF (Cf==3) THEN !SAMPSON
    
    F1= -qaux1*C(i,j)/(1.0D0+dt*C(i,j))
    F2= -qaux2*C(i,j)/(1.0D0+dt*C(i,j))

    ELSE IF (Cf==4) THEN ! CONSTANT Cf
    
    F1= -alpha1*C(i,j)*Fr2/(1.0D0+2.0D0*dt*C(i,j)*sqrt(qaux1**2.0D0+qaux2**2.0D0)/(q(1,i,j)**2.0D0))
    F2= -alpha2*C(i,j)*Fr2/(1.0D0+2.0D0*dt*C(i,j)*sqrt(qaux1**2.0D0+qaux2**2.0D0)/(q(1,i,j)**2.0D0))
    
    END IF


END IF

!Limitador para F

Flim1=-1.0D0*qaux1/dt
Flim2=-1.0D0*qaux2/dt
IF (qaux1>=0.0D0.AND.F1<Flim1) then
F1=Flim1
END IF

IF (qaux2>=0.0D0.AND.F2<Flim2) then
F2=Flim2
END IF

IF (qaux1<=0.0D0.AND.F1>Flim1) then
F1=Flim1
END IF

IF (qaux2<=0.0D0.AND.F2>Flim2) then
F2=Flim2
END IF



qfaux(1,i,j)=q(1,i,j)
qfaux(2,i,j)=qaux1+dt*F1
qfaux(3,i,j)=qaux2+dt*F2



!Tolerancia
  IF (qfaux(1,i,j).le.hmin) THEN	!Dry Cell
    qfaux(1,i,j)=0.0D0
    qfaux(2,i,j)=0.0D0
    qfaux(3,i,j)=0.0D0
  END IF
  
  qf(1,i,j)=qfaux(1,i,j)
  
  IF (qf(1,i,j)==0.0D0) THEN
  qf(2,i,j)=qfaux(2,i,j)
  qf(3,i,j)=qfaux(3,i,j)
  
  ELSE
  qf(2,i,j)=qfaux(2,i,j)/qfaux(1,i,j)
  qf(3,i,j)=qfaux(3,i,j)/qfaux(1,i,j)
  END IF
  
  IF (abs(qf(2,i,j)).le.kappa) THEN
  qf(2,i,j)=0.0D0
  END IF
  
  IF (abs(qf(3,i,j)).le.kappa) THEN
  qf(3,i,j)=0.0D0
  END IF

END DO; END DO

END SUBROUTINE friccion
