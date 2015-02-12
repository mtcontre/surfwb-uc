SUBROUTINE VyC
!Subroutine for the CFL condition
!o adimensional
!
USE global_variables
implicit none


integer	:: i,j
real (kind=8)	::  U1, U2
! real (kind=8), dimension(2,Nbx,Nby)	:: qnew_global_abs
!ya calculado todo, 
DO i=1,Nbx; DO j=1,Nby
! 	qnew_global_abs(1,i,j)=abs(qnew_global(2,i,j))
! 	qnew_global_abs(2,i,j)=abs(qnew_global(3,i,j))
	V_global(i,j)=sqrt(qnew_global(2,i,j)**2.0D0+qnew_global(3,i,j)**2.0D0)!???se usa este ???	
	U1=qnew_global(2,i,j)*xi_global(1,i,j)+qnew_global(3,i,j)*xi_global(2,i,j)
	U2=qnew_global(2,i,j)*eta_global(1,i,j)+qnew_global(3,i,j)*eta_global(2,i,j)
	C_global(i,j)=sqrt(qnew_global(1,i,j)/FR2)
	VC(i,j)=V_global(i,j)+C_global(i,j)
	S1_global(i,j)=abs(U1)+C_global(i,j)*sqrt(xi_global(1,i,j)**2+xi_global(2,i,j)**2)
	S2_global(i,j)=abs(U2)+C_global(i,j)*sqrt(eta_global(1,i,j)**2+eta_global(2,i,j)**2)
	
END DO; END do
END SUBROUTINE VyC

subroutine checkSbc(Sg1u1,Sg1u2,Sg2u1,Sg2u2,varsg1,varsg2,metrg1,metrg2)
!Calculates Riemann invariant propagation speed for boundary conditions g1,g2
!to calculate dt according to numerical stability considering boundary values
  use global_variables
  real (kind=8) 		:: Sg1u1,Sg1u2,Sg2u1,Sg2u2 !output
  real (kind=8),dimension(3) :: varsg1,varsg2			!input hydr. vars	
  real (kind=8),dimension(4) :: metrg1,metrg2			!input metrics for bc
  real (kind=8)			:: hg1,ug1,vg1,hg2,ug2,vg2	!hydr vars to get from 'varsg1/2'
  real (kind=8)			:: xig1x,xig1y,etag1x,etag1y	!metrics from metrg1
  real (kind=8)			:: xig2x,xig2y,etag2x,etag2y	!metrics from metrg2
  real (kind=8)			:: U1,U2,Cg			!Contravariant speeds, adim. celerity  
  hg1=varsg1(1)
  ug1=varsg1(2)
  vg1=varsg1(3)
  hg2=varsg2(1)
  ug2=varsg2(2)
  vg2=varsg2(3)
  xig1x=metrg1(1)
  xig1y=metrg1(2)
  etag1x=metrg1(3)
  etag1y=metrg1(4)
  xig2x=metrg1(1)
  xig2y=metrg1(2)
  etag2x=metrg1(3)
  etag2y=metrg1(4)
  
  U1=ug1*xig1x +vg1*xig1y 
  U2=ug1*etag1x+vg1*etag1y
  Cg=sqrt(hg1/FR2)
  Sg1u1=abs(U1)+Cg*sqrt(xig1x**2+xig2x**2)
  Sg1u2=abs(U2)+Cg*sqrt(etag1x**2+etag1y**2)
	  
  U1=ug2*xig1x+vg2*xig1y 
  U2=ug2*etag1x+vg2*etag1y
  Cg=sqrt(hg2/FR2)
  Sg2u1=abs(U1)+Cg*sqrt(xig1x**2+xig1y**2)
  Sg2u2=abs(U2)+Cg*sqrt(etag1x**2+etag1y**2)	  
end subroutine