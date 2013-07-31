SUBROUTINE solver2
!2nd order Runge-Kutta time integration
! NO friction

USE global_variables
USE geometries
implicit none

!Local variables definition
integer	:: i,j
real (kind=8)	:: u1max, v1max, u2max, v2max, VCmax, Vmax,maxC, umax, vvmax, dxi1, deta1, kappa2, maxDF, maxDFloc
real (kind=8),dimension(2):: locu1, locv1, DFloc
real (kind=8),dimension(:,:,:),allocatable :: q1,q2, &
				     q1T,q2T, &
				     xi_T, eta_T,zT, &
				     F1mas,F2mas,&
				     G1mas,G2mas,&
				     F1menos,F2menos, &
				     G1menos,G2menos, &
				     qaux, qaux1, qaux2, &
				     SC, k1, k2					!T Includes ghost cells
				     !F1,F2,F3,F4, &
				     !G1,G2,G3,G4, &	

real (kind=8),dimension(:), allocatable	:: deltaF, dh
allocate(deltaF(Nbx),dh(Nbx))

!Allocate intermediate state variables
allocate( q1(3,Nbx,Nby),q2(3,Nbx,Nby), &
	  F1mas(3,Nbx,Nby),F2mas(3,Nbx,Nby), &
	  G1mas(3,Nbx,Nby),G2mas(3,Nbx,Nby),&
	  F1menos(3,Nbx,Nby),F2menos(3,Nbx,Nby), &
	  G1menos(3,Nbx,Nby),G2menos(3,Nbx,Nby), &
	  qaux(3,Nbx,Nby),qaux1(3,Nbx,Nby),qaux2(3,Nbx,Nby), &
	  SC(3,Nbx,Nby), k1(3,Nbx,Nby), k2(3,Nbx,Nby))
	  
	  !F1(3,Nbx,Nby),F2(3,Nbx,Nby),F3(3,Nbx,Nby),F4(3,Nbx,Nby), &
	  !G1(3,Nbx,Nby),G2(3,Nbx,Nby),G3(3,Nbx,Nby),G4(3,Nbx,Nby), &

!Allocate intermediate state variables with ghost cells
allocate(q1T(3,Nbx+4,Nby+4),q2T(3,Nbx+4,Nby+4))
allocate(xi_T(2,Nbx+4,Nby+4),eta_T(2,Nbx+4,Nby+4),zT(1,Nbx+4,Nby+4))

kappa2=1e-10
!q=(h,u,v)



!I need qaux=(h,hu,hv) for the time integration, and q1 must have (h,u,v) for the bcs and fluxes functions:

DO i=1,Nbx; DO j=1,Nby
  qaux(1,i,j)=qold_global(1,i,j)
  qaux(2,i,j)=qold_global(1,i,j)*qold_global(2,i,j)
  qaux(3,i,j)=qold_global(1,i,j)*qold_global(3,i,j)
END DO; END DO

!Cleaning qaux

Do i=1,Nbx; Do j=1,Nby
  if (qaux(1,i,j).le.hmin) then
  qaux(1,i,j)=0.0D0
  qaux(2,i,j)=0.0D0
  qaux(3,i,j)=0.0D0
  end if
  
  if (abs(qaux(2,i,j)).le.kappa) then
  qaux(2,i,j)=0.0D0
  end if
  
  if (abs(qaux(3,i,j)).le.kappa) then
  qaux(3,i,j)=0.0D0
  end if
  
End do; End do

!First RK Stage

!Calculates q(n+1/2*)

call bcs(fopt,Cf,MCoef,1,t,0.5D0*dt,FR2,caso,qold_global,z_global,Nbx,Nby,CB,xi_global,eta_global,aj_global,dxi,deta,q1T,xi_T,eta_T,zT)

call fluxes(CB,mmopt,hmin,q1T,zT,xi_T,eta_T,dxi,deta,Nbx,Nby,FR2,F1mas,F1menos,G1mas,G1menos,SC)

DO i=1,Nbx; DO j=1,Nby
k1(:,i,j)=(SC(:,i,j)-(aj_global(i,j)/dxi)*(F1menos(:,i,j)-F1mas(:,i,j))-(aj_global(i,j)/deta)*(G1menos(:,i,j)-G1mas(:,i,j)))
qaux1(:,i,j)=qaux(:,i,j)+0.5D0*dt*k1(:,i,j)

END DO; END DO

!q1
DO i=1,Nbx; DO j=1,Nby

  IF (qaux1(1,i,j).le.hmin) THEN	!Dry Cell
    qaux1(1,i,j)=0.0D0
    qaux1(2,i,j)=0.0D0
    qaux1(3,i,j)=0.0D0
  END IF
  
  q1(1,i,j)=qaux1(1,i,j)
  
  IF (q1(1,i,j)==0.0D0) THEN
  q1(2,i,j)=qaux1(2,i,j)
  q1(3,i,j)=qaux1(3,i,j)
  
  ELSE
  q1(2,i,j)=qaux1(2,i,j)/qaux1(1,i,j)
  q1(3,i,j)=qaux1(3,i,j)/qaux1(1,i,j)
  END IF
  
  IF (abs(q1(2,i,j)).le.kappa) THEN
  q1(2,i,j)=0.0D0
  END IF
  
  IF (abs(q1(3,i,j)).le.kappa) THEN
  q1(3,i,j)=0.0D0
  END IF
  
END DO; END DO


!Second RK Stage
!Calculates q(n+1)

call bcs(fopt,Cf,MCoef,2,t+0.5D0*dt,0.5D0*dt,FR2,caso,q1,z_global,Nbx,Nby,CB,xi_global,eta_global,aj_global,dxi,deta,q2T,xi_T,eta_T,zT)

call fluxes(CB,mmopt,hmin,q2T,zT,xi_T,eta_T,dxi,deta,Nbx,Nby,FR2,F2mas,F2menos,G2mas,G2menos,SC)

DO i=1,Nbx; DO j=1,Nby
k2(:,i,j)=(SC(:,i,j)-aj_global(i,j)/dxi*(F2menos(:,i,j)-F2mas(:,i,j))-aj_global(i,j)/deta*(G2menos(:,i,j)-G2mas(:,i,j)))

qaux2(:,i,j)=qaux(:,i,j)+dt*k2(:,i,j)
END DO; END DO


!q2
DO i=1,Nbx; DO j=1,Nby

  IF (qaux2(1,i,j)<=hmin) THEN	!Dry Cell
    qaux2(1,i,j)=0.0D0
    qaux2(2,i,j)=0.0D0
    qaux2(3,i,j)=0.0D0
   
  END IF
  
  q2(1,i,j)=qaux2(1,i,j)
  
  IF (q2(1,i,j)==0.0D0) THEN
  q2(2,i,j)=0.0D0
  q2(3,i,j)=0.0D0
  
  ELSE
  q2(2,i,j)=qaux2(2,i,j)/qaux2(1,i,j)
  q2(3,i,j)=qaux2(3,i,j)/qaux2(1,i,j)
  END IF
  
  IF (abs(q2(2,i,j)).le.kappa) THEN
  q2(2,i,j)=0.0D0
  END IF
  
  IF (abs(q2(3,i,j)).le.kappa) THEN
  q2(3,i,j)=0.0D0
  END IF
  
END DO; END DO





!New Result (n+1)
qnew_global=q2

!pause
!V y C para calcular dt
call VyC2

!-----------------------------------------

!-----------------------------------------
END SUBROUTINE solver2



SUBROUTINE VyC2
!Subroutine for the CFL condition
!o adimensional
USE global_variables
implicit none


integer	:: i,j
integer, dimension(2)::locMaxV
real (kind=8)	:: raiz, u2, v2, maxC, maxV, maxu1, maxu2, U1, Um2
real (kind=8), dimension(2,Nbx,Nby)	:: qnew_global_abs
DO i=1,Nbx; DO j=1,Nby
	u2=qnew_global(2,i,j)**2.0D0
	v2=qnew_global(3,i,j)**2.0D0
	qnew_global_abs(1,i,j)=abs(qnew_global(2,i,j))
	qnew_global_abs(2,i,j)=abs(qnew_global(3,i,j))
	raiz=u2+v2
	V_global(i,j)=sqrt(raiz)
	C_global(i,j)=sqrt(qnew_global(1,i,j)/FR2)
	VC(i,j)=V_global(i,j)+C_global(i,j)
	
	U1=qold_global(2,i,j)*xi_global(1,i,j)+qold_global(3,i,j)*xi_global(2,i,j)
	Um2=qold_global(2,i,j)*eta_global(1,i,j)+qold_global(3,i,j)*eta_global(2,i,j)

	S1_global(i,j)=U1+C_global(i,j)*sqrt(xi_global(1,i,j)**2+xi_global(2,i,j)**2)
	S2_global(i,j)=Um2+C_global(i,j)*sqrt(eta_global(1,i,j)**2+eta_global(2,i,j)**2)
	
END DO; END DO


END SUBROUTINE VyC2
