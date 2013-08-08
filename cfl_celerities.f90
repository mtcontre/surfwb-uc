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


! subroutine VyC!stability_celerities_in
!   !allocates celerities for in-domain cells to  calculate dt (later in main.f90)
!   !accordin to cfl stability condition
! 
!   use global_variables
!   use geometries
!   use coords
!   use custombc
!   integer ::i,j
!   real (kind=8) :: U1, U2,Cg
!   !Velocidades Contravariantes para CFL
!   !JGM: i guess max. riemann invariant propagation speed shoul be |u|+c, and not |u+c| nor u+c
!   !modified 27/6/2013
!   do i=1,Nbx; do j=1,Nby
! 	  U1=qold_global(2,i,j)*xi_global(1,i,j)+qold_global(3,i,j)*xi_global(2,i,j)
! 	  U2=qold_global(2,i,j)*eta_global(1,i,j)+qold_global(3,i,j)*eta_global(2,i,j)
! 	  S1_global(i,j)=abs(U1)+C_global(i,j)*sqrt(xi_global(1,i,j)**2+xi_global(2,i,j)**2)
! 	  S2_global(i,j)=abs(U2)+C_global(i,j)*sqrt(eta_global(1,i,j)**2+eta_global(2,i,j)**2)
!   end do; end do
!   
!   
! end subroutine


subroutine stability_celerities_boundary_init!creo 2 para no tener conflicto con las variables, este es solo para el instante inicial
  !allocates celerities for in-boundary cells to calculate dt (later in main.f90) 
  !accordin to cfl stability condition
  
  use global_variables
  use geometries
  use coords
  use custombc
  integer ::i,j
  real (kind=8) :: U1, U2,Cg

  if (flagxi0.eq.1) then !if xi0 bound ==custom    
    allocate(Sxi0(Nby,4))
    do j=1,Nby
!     checkSbc(Sg1u1,Sg1u2,Sg2u1,Sg2u2,varsg1,varsg2,metrg1,metrg2)
      call checkSbc(   Sxi0(j,1),Sxi0(j,2),Sxi0(j,3),Sxi0(j,4)    ,&
	   (/qxi0g1(2,j,1),qxi0g1(3,j,1),qxi0g1(4,j,1)/), &
	   (/qxi0g2(2,j,1),qxi0g2(3,j,1),qxi0g2(4,j,1)/), &
	   (/xi_global(1,2,j),xi_global(2,2,j),eta_global(1,2,j),eta_global(2,2,j)/), &
	   (/xi_global(1,1,j),xi_global(2,1,j),eta_global(1,1,j),eta_global(2,1,j)/))
      !a la ghost cell 1 le corresponde la metrica de la celda 2 (las metricas se copian simetricas)
    end do
 end if
  
  if (flagxiN.eq.1) then !if xi0 bound ==custom
    allocate(SxiN(Nby,4))
    do j=1,Nby
      call checkSbc(   SxiN(j,1),SxiN(j,2),SxiN(j,3),SxiN(j,4)    ,&
	   (/qxiNg1(2,j,1),qxiNg1(3,j,1),qxiNg1(4,j,1)/), &
	   (/qxiNg2(2,j,1),qxiNg2(3,j,1),qxiNg2(4,j,1)/), &
	   (/xi_global(1,Nbx-1,j),xi_global(2,Nbx-1,j),eta_global(1,Nbx-1,j),eta_global(2,Nbx-1,j)/), &
	   (/xi_global(1,Nbx,j),xi_global(2,Nbx,j),eta_global(1,Nbx,j),eta_global(2,Nbx,j)/))
    end do
  end if
  
  if (flageta0.eq.1) then
    allocate(Seta0(Nbx,4))
    do i=1,Nbx
      call checkSbc(   Seta0(j,1),Seta0(j,2),Seta0(j,3),Seta0(j,4)    ,&
	   (/qeta0g1(2,j,1),qeta0g1(3,j,1),qeta0g1(4,j,1)/), &
	   (/qeta0g2(2,j,1),qeta0g2(3,j,1),qeta0g2(4,j,1)/), &
	   (/xi_global(1,i,2),xi_global(2,i,2),eta_global(1,i,2),eta_global(2,i,2)/), &
	   (/xi_global(1,i,1),xi_global(2,i,1),eta_global(1,i,1),eta_global(2,i,1)/))
    end do
  end if

  if (flagetaN.eq.1) then
    allocate(SetaN(Nbx,4))
    do i=1,Nbx
      call checkSbc(   SetaN(j,1),SetaN(j,2),SetaN(j,3),SetaN(j,4)    ,&
	   (/qetaNg1(2,j,1),qetaNg1(3,j,1),qetaNg1(4,j,1)/), &
	   (/qetaNg2(2,j,1),qetaNg2(3,j,1),qetaNg2(4,j,1)/), &
	   (/xi_global(1,i,Nby),xi_global(2,i,Nby),eta_global(1,i,Nby),eta_global(2,i,Nby)/), &
	   (/xi_global(1,i,Nby-1),xi_global(2,i,Nby-1),eta_global(1,i,Nby-1),eta_global(2,i,Nby-1)/))
    end do
  end if


end subroutine

subroutine stability_celerities_boundary(qt) !este es para instantes distintos al inicial
  !allocates celerities for in-boundary cells to calculate dt (later in main.f90) 
  !accordin to cfl stability condition
  
  use global_variables
  use geometries
  use coords
  use custombc
  integer ::i,j
  real (kind=8) :: U1, U2,Cg
  real (kind=8), dimension(3,Nbx+4,Nby+4)	::qt
  if (flagxi0.eq.1) then !if xi0 bound ==custom    
    do j=1,Nby
!     checkSbc(Sg1u1,Sg1u2,Sg2u1,Sg2u2,varsg1,varsg2,metrg1,metrg2)
      call checkSbc(   Sxi0(j,1),Sxi0(j,2),Sxi0(j,3),Sxi0(j,4)    ,&
	   (/qt(1,1,j+2),qt(2,1,j+2),qt(3,1,j+2)/), &
	   (/qt(1,2,j+2),qt(2,2,j+2),qt(3,2,j+2)/), &
	   (/xi_global(1,2,j),xi_global(2,2,j),eta_global(1,2,j),eta_global(2,2,j)/), &
	   (/xi_global(1,1,j),xi_global(2,1,j),eta_global(1,1,j),eta_global(2,1,j)/))
      !a la ghost cell 1 le corresponde la metrica de la celda 2 (las metricas se copian simetricas)
    end do
 end if
  
  if (flagxiN.eq.1) then !if xi0 bound ==custom
    do j=1,Nby
      call checkSbc(   SxiN(j,1),SxiN(j,2),SxiN(j,3),SxiN(j,4)    ,&
	   (/qt(1,Nx+3,j+2),qt(2,Nx+3,j+2),qt(3,Nx+3,j+2)/), &
	   (/qt(1,Nx+4,j+2),qt(2,Nx+4,j+2),qt(3,Nx+4,j+2)/), &
	   (/xi_global(1,Nbx-1,j),xi_global(2,Nbx-1,j),eta_global(1,Nbx-1,j),eta_global(2,Nbx-1,j)/), &
	   (/xi_global(1,Nbx,j),xi_global(2,Nbx,j),eta_global(1,Nbx,j),eta_global(2,Nbx,j)/))
    end do
  end if
  
  if (flageta0.eq.1) then
    do i=1,Nbx
      call checkSbc(   Seta0(j,1),Seta0(j,2),Seta0(j,3),Seta0(j,4)    ,&
	   (/qt(1,i+2,1),qt(2,i+2,1),qt(3,i+2,1)/), &
	   (/qt(1,i+2,2),qt(2,i+2,2),qt(3,i+2,2)/), &
	   (/xi_global(1,i,2),xi_global(2,i,2),eta_global(1,i,2),eta_global(2,i,2)/), &
	   (/xi_global(1,i,1),xi_global(2,i,1),eta_global(1,i,1),eta_global(2,i,1)/))
    end do
  end if

  if (flagetaN.eq.1) then
    do i=1,Nbx
      call checkSbc(   SetaN(j,1),SetaN(j,2),SetaN(j,3),SetaN(j,4)    ,&
	   (/qt(1,i+2,Ny+3),qt(2,i+2,Ny+3),qt(3,i+2,Ny+3)/), &
	   (/qt(1,i+2,Ny+4),qt(2,i+2,Ny+4),qt(3,i+2,Ny+4)/), &
	   (/xi_global(1,i,Nby),xi_global(2,i,Nby),eta_global(1,i,Nby),eta_global(2,i,Nby)/), &
	   (/xi_global(1,i,Nby-1),xi_global(2,i,Nby-1),eta_global(1,i,Nby-1),eta_global(2,i,Nby-1)/))
    end do
  end if


end subroutine
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