SUBROUTINE FLUXES(CB,mmopt,hmin,qt,zt,xit,etat,dxi,deta,Nx,Ny,FR2,Fmas,Fmenos,Gmas,Gmenos,SC)
!Funcion que calcula los flujos numericos a traves de las interfaces de cada celda
!OFICIAL!!! ESTE ES EL QUE SE USA

use senales
use Jacobianos
implicit none

integer :: Nx,Ny, mmopt

real (kind=8), dimension(3,Nx+4,Ny+4)	::qt 
real (kind=8), dimension(Nx+4,Ny+4)	::zt, hzt
real (kind=8), dimension(2,Nx+4,Ny+4)	::xit
real (kind=8), dimension(2,Nx+4,Ny+4)	::etat
real (kind=8), dimension(2,Nx+2,Ny)	::xi
real (kind=8), dimension(2,Nx,Ny+2)	::eta
real (kind=8), dimension(Nx+1,Ny)	::Jac1
real (kind=8), dimension(Nx,Ny+1)	::Jac2

real (kind=8), dimension(2,Nx+1,Ny)	::xim
real (kind=8), dimension(2,Nx,Ny+1)	::etam

real (kind=8), dimension(3,Nx,Ny)	::q1, SC, SF, SG
real (kind=8), dimension(3,Nx,Ny)	::Fmas,Fmenos,Gmas,Gmenos 
real (kind=8), dimension(3,Nx+1,Ny)	::QiL, QiR, Qstari
real (kind=8), dimension(3,Nx,Ny+1)	::QjL, QjR, Qstarj
real (kind=8), dimension(Nx+1,Ny)	::HziL,HziR,&
				  ZiL,ZiR, &
				  hmasi,hmenosi
real (kind=8), dimension(Nx,Ny+1)	::HzjL,HzjR,&
				  ZjL,ZjR, &
				  hmasj,hmenosj	
integer::i,j,k
real (kind=8),dimension(3) :: x,y,z,mm, Fs,Gs, qmas, qmenos, Fb,Gb, qs, QauxiL, QauxiR,maxm1, maxm2
real (kind=8),dimension(3,3) :: Qauxi, Qauxj
real (kind=8) :: FR2,x1,y1,z1,mm1, zm, zmj, L, H, hmin, dxi, deta, hzmin, U1, U2
real (kind=8),dimension(2) ::hzminloc, ximas, ximenos, etamas, etamenos

!Senales entrando al dominio
integer,dimension(4)::CB

!1.Minmod reconstruction of the state variables

!Reconstruction of h,u y v (not h, hu, hv because ur,ul,vr,vl are needed)
!For every xi(i) i just need real (kind=8) etas without ghost cells, j between (3 y Ny+2) of qt matrix
!For every eta(i) i just need real (kind=8) xis without ghost cells, i between (3 y Nx+2) of qt matrix

!Estoy calculando desde el i=1 y necesito i=0R tb 
!Estoy calculando hasta el i=Nx y necesito tb i=Nx+1
!Osea desde i=2 a i=Nx+3


!$omp parallel private(x, y, z, mm, maxm1, maxm2, &
!$omp x1, y1, z1, mm1,  zm, zmj, &
!$omp qmenos, ximenos, qmas, ximas, qs, etamenos, etamas, &
!$omp U1, Fs, Fb, U2, Gs, Gb)

!---------loop 1---------------x,y,z,mm,maxm1,maxm2
!$omp do
Do i=1,Nx+1; Do j=1,Ny
  !Se necesita QiL de la celda 1 a la N+1
  !Para calcularlo necesito tb i=0 y N+2
  !En términos de qt se usara i=2 a N+3
  x(:)=(qt(:,i+2,j+2)-qt(:,i+1,j+2))*(1.0D0/dxi)	!Q1-Q0=Qt3-Qt2		!QN+1-QN=QtN+3-QtN+2	!Qi-Qi-1
  y(:)=(qt(:,i+3,j+2)-qt(:,i+2,j+2))*(1.0D0/dxi)	!Q2-Q1=Qt4-Qt3		!QN+2-QN-1=QtN+4-QtN+3	!Qi+1-Qi
  z(:)=(qt(:,i+3,j+2)-qt(:,i+1,j+2))/(2.0D0*dxi)	!Qi+1-Qi-1


  Do k=1,3
    if (mmopt==1) then
      call minmod(x(k),y(k),mm(k))
    else if (mmopt==2) then
      call minmod((2.0D0*x(k)),y(k),maxm1(k))
      call minmod(x(k),(2.0D0*y(k)),maxm2(k))
      call maxmod(maxm1(k),maxm2(k),mm(k))
    else if (mmopt==3) then ! Double minmod
      call DOUBLEMINMOD(x(k),y(k),z(k),mm(k))
    else
      call mc(2.0D0*x(k),2.0D0*y(k),z(k),mm(k))      
    end if
  end Do
  
  QiL(:,i,j)=qt(:,i+2,j+2)-0.5D0*dxi*mm(:)	!Q1L=Q1-dxi/2*minmod(x,y)
  
  
  
  !Se necesita QiR de la celda 0 a N
  !Para calcularlo necesito usar i=-1 y N+1
  !En términos de qt se usará: i=1, a N+3
  
  x(:)=(qt(:,i+1,j+2)-qt(:,i,j+2))*(1.0D0/dxi)	!Q0-Q-1=Qt2-Qt1		!QN-QN-1=QtN+2-QtN+1
  y(:)=(qt(:,i+2,j+2)-qt(:,i+1,j+2))*(1.0D0/dxi)	!Q1-Q0=Qt3-Qt2		!QN+1-QN=QtN+3-QtN+2
  z(:)=(qt(:,i+2,j+2)-qt(:,i,j+2))*(1.0D0/(2.0D0*dxi))	
  
  Do k=1,3
    if (mmopt==1) then
      call minmod(x(k),y(k),mm(k))
    else if (mmopt==2) then
      call minmod((2.0D0*x(k)),y(k),maxm1(k))
      call minmod(x(k),(2.0D0*y(k)),maxm2(k))
      call maxmod(maxm1(k),maxm2(k),mm(k))
    else if (mmopt==3) then ! Double minmod
      call DOUBLEMINMOD(x(k),y(k),z(k),mm(k))
    else
      call mc(2.0D0*x(k),2.0D0*y(k),z(k),mm(k))    
    end if

  end Do

  QiR(:,i,j)=qt(:,i+1,j+2)+0.5D0*dxi*mm(:)
  
  !Reemplazo valores del borde si es que hay bordes gen/abs
  !Igual calcula antes un valor falso usando todas las celdas ficticias que tienen valores, -1=qA, 0=qA
  IF (i==1.AND.(CB(1)==4.OR.CB(1)==5)) THEN
  QiR(1,1,j)=qA1(1,j)
  QiR(2,1,j)=qA1(2,j)
  QiR(3,1,j)=qA1(3,j)
  !print*,'QiR=',QiR(1,1,j), QiR(2,1,j), QiR(3,1,j), j
  !print*,'qA1=',qA1(1,j), qA1(2,j), qA1(3,j), j
  END IF

  IF (i==Nx+1.AND.(CB(2)==4.OR.CB(2)==5)) THEN
  !print*, 'holo'
  QiL(1,Nx+1,j)=qA2(1,j)
  QiL(2,Nx+1,j)=qA2(2,j)
  QiL(3,Nx+1,j)=qA2(3,j)

  END IF
  
 
End Do; End Do
!$omp end do NOWAIT

!---------loop 2---------------x,y,z,mm,maxm1,maxm2
!$omp do
Do i=1,Nx; Do j=1,Ny+1
  
  !Se necesita QjL de la celda 1 a la N+1
  !Para calcularlo necesito tb j=0 y N+2
  !En términos de qt se usara j=2 a N+3
  
  x(:)=(qt(:,i+2,j+2)-qt(:,i+2,j+1))/deta
  y(:)=(qt(:,i+2,j+3)-qt(:,i+2,j+2))/deta
  z(:)=(qt(:,i+2,j+3)-qt(:,i+2,j+1))/(2.0D0*deta)
  Do k=1,3
  if (mmopt==1) then
    call minmod(x(k),y(k),mm(k))
  else if (mmopt==2) then
    call minmod((2.0D0*x(k)),y(k),maxm1(k))
    call minmod(x(k),(2.0D0*y(k)),maxm2(k))
    call maxmod(maxm1(k),maxm2(k),mm(k))
  else if (mmopt==3) then ! Double minmod
      call DOUBLEMINMOD(x(k),y(k),z(k),mm(k))
  else
    call mc(2.0D0*x(k),2.0D0*y(k),z(k),mm(k))    
  end if
  end Do
  
  QjL(:,i,j)=qt(:,i+2,j+2)-0.5D0*deta*mm(:)
  
  !Se necesita QjR de la celda 0 a N
  !Para calcularlo necesito usar j=-1 y N+1
  !En términos de qt se usará: j=1, a N+3

  x(:)=(qt(:,i+2,j+1)-qt(:,i+2,j))/deta
  y(:)=(qt(:,i+2,j+2)-qt(:,i+2,j+1))/deta
  z(:)=(qt(:,i+2,j+2)-qt(:,i+2,j))/(2.0D0*deta)
  
  Do k=1,3
  if (mmopt==1) then
    call minmod(x(k),y(k),mm(k))
  else if (mmopt==2) then
    call minmod((2.0D0*x(k)),y(k),maxm1(k))
    call minmod(x(k),(2.0D0*y(k)),maxm2(k))
    call maxmod(maxm1(k),maxm2(k),mm(k))
  else if (mmopt==3) then ! Double minmod
      call DOUBLEMINMOD(x(k),y(k),z(k),mm(k))
  else
    call mc(2.0D0*x(k),2.0D0*y(k),z(k),mm(k))    
  end if
  
  end Do
  
  QjR(:,i,j)=qt(:,i+2,j+1)+0.5D0*deta*mm(:)
  
  !Reemplazo valores del borde si es que hay bordes gen/abs
  !Igual calcula antes un valor falso usando todas las celdas ficticias que tienen valores, -1=0, 0=qA
  IF (j==1.AND.(CB(3)==4.OR.CB(3)==5)) THEN
  QjR(1,i,1)=qA3(1,i)
  QjR(2,i,1)=qA3(2,i)
  QjR(3,i,1)=qA3(3,i)
  END IF
  IF (j==Ny+1.AND.(CB(4)==4.OR.CB(4)==5)) THEN
  QjL(1,i,Ny+1)=qA4(1,i)
  QjL(2,i,Ny+1)=qA4(2,i)
  QjL(3,i,Ny+1)=qA4(3,i)
  END IF
  
End Do; End Do
!$omp end do NOWAIT

!2.Linear Reconstruction of the free surface

!---------loop 3--------------- no privates
!$omp do
DO i=1,(Nx+4); DO j=1,(Ny+4)
  hzt(i,j)=qt(1,i,j)+zt(i,j)
END DO; END DO
!$omp end do

!---------loop 4---------------x1,y1,z1,mm1,maxm1,maxm2
!$omp do
Do i=1,Nx+1; Do j=1,Ny

  !Se necesita HZiL de la celda 1 a la N+1
  !Para calcularlo necesito tb i=0 y N+2
  !En términos de qt se usara i=2 a N+3
    
  x1=(hzt(i+2,j+2)-hzt(i+1,j+2))/dxi
  y1=(hzt(i+3,j+2)-hzt(i+2,j+2))/dxi
  z1=(hzt(i+3,j+2)-hzt(i+1,j+2))/(2.0D0*dxi)
  if (mmopt==1) then
    call minmod(x1,y1,mm1)
  else if (mmopt==2) then
    call minmod((2.0D0*x1),y1,maxm1(1))
    call minmod(x1,(2.0D0*y1),maxm2(1))
    call maxmod(maxm1(1),maxm2(1),mm1)
  else if (mmopt==3) then ! Double minmod
      call DOUBLEMINMOD(x1,y1,z1,mm1)
  else
    call mc(2.0D0*x1,2.0D0*y1,z1,mm1)
  end if
  
  HziL(i,j)=hzt(i+2,j+2)-0.5D0*dxi*mm1
  
  
  
  !Se necesita HZiR de la celda 0 a N
  !Para calcularlo necesito usar i=-1 y N+1
  !En términos de qt se usará: i=1, a N+3
  
  x1=(hzt(i+1,j+2)-hzt(i,j+2))/dxi
  y1=(hzt(i+2,j+2)-hzt(i+1,j+2))/dxi
  z1=(hzt(i+2,j+2)-hzt(i,j+2))/(2.0D0*dxi)
  
  if (mmopt==1) then
    call minmod(x1,y1,mm1)
  else if (mmopt==2) then
    call minmod((2.0D0*x1),y1,maxm1(1))
    call minmod(x1,(2.0D0*y1),maxm2(1))
    call maxmod(maxm1(1),maxm2(1),mm1)
  else if (mmopt==3) then ! Double minmod
      call DOUBLEMINMOD(x1,y1,z1,mm1)
  else
    call mc(2.0D0*x1,2.0D0*y1,z1,mm1)
  end if
  
  HziR(i,j)=hzt(i+1,j+2)+0.5D0*dxi*mm1

  !Reemplazo valores del borde si es que hay bordes gen/abs
  !Igual calcula antes un valor falso usando todas las celdas ficticias que tienen valores, -1=0, 0=qA
  IF (i==1.AND.(CB(1)==4.OR.CB(1)==5)) THEN
  HziR(1,j)=qA1(1,j)+zA1(j)
  END IF
  IF (i==Nx+1.AND.(CB(2)==4.OR.CB(2)==5)) THEN
  HziL(Nx+1,j)=qA2(1,j)+zA2(j)
  END IF
  
  
End Do; End Do
!$omp end do NOWAIT

!---------loop 5---------------x1,y1,z1,mm1,maxm1,maxm2
!$omp do
Do i=1,Nx; Do j=1,Ny+1

  !Se necesita HZjL de la celda 1 a la N+1
  !Para calcularlo necesito tb j=0 y N+2
  !En términos de qt se usara j=2 a N+3
  
  x1=(hzt(i+2,j+2)-hzt(i+2,j+1))/deta
  y1=(hzt(i+2,j+3)-hzt(i+2,j+2))/deta
  z1=(hzt(i+2,j+3)-hzt(i+2,j+1))/(2.0D0*deta)
  select case(mmopt)
     case(1)
      call minmod(x1,y1,mm1)
     case(2)
      call minmod((2.0D0*x1),y1,maxm1(1))
      call minmod(x1,(2.0D0*y1),maxm2(1))
      call maxmod(maxm1(1),maxm2(1),mm1)
     case(3)
      call DOUBLEMINMOD(x1,y1,z1,mm1)
     case default
      call mc(2.0D0*x1,2.0D0*y1,z1,mm1)
  end select
  
!   if (mmopt==1) then
!     call minmod(x1,y1,mm1)
!   else if (mmopt==2) then
!     call minmod((2.0D0*x1),y1,maxm1(1))
!     call minmod(x1,(2.0D0*y1),maxm2(1))
!     call maxmod(maxm1(1),maxm2(1),mm1)
!   else if (mmopt==3) then ! Double minmod
!       call DOUBLEMINMOD(x1,y1,z1,mm1)
!   else
!     call mc(2.0D0*x1,2.0D0*y1,z1,mm1)
!   end if
  
  HzjL(i,j)=hzt(i+2,j+2)-0.5D0*deta*mm1
  
  !Se necesita HZjR de la celda 0 a N
  !Para calcularlo necesito usar j=-1 y N+1
  !En términos de qt se usará: j=1, a N+3

  x1=(hzt(i+2,j+1)-hzt(i+2,j))/deta
  y1=(hzt(i+2,j+2)-hzt(i+2,j+1))/deta
  z1=(hzt(i+2,j+2)-hzt(i+2,j))/(2.0D0*deta)
  
  select case (mmopt)
    case(1)
      call minmod(x1,y1,mm1)
    case(2)
      call minmod((2.0D0*x1),y1,maxm1(1))
      call minmod(x1,(2.0D0*y1),maxm2(1))
      call maxmod(maxm1(1),maxm2(1),mm1)
    case(3)
      call DOUBLEMINMOD(x1,y1,z1,mm1)
    case default
      call mc(2.0D0*x1,2.0D0*y1,z1,mm1)
  end select
!   if (mmopt==1) then
!     call minmod(x1,y1,mm1)
!   else if (mmopt==2) then
!     call minmod((2.0D0*x1),y1,maxm1(1))
!     call minmod(x1,(2.0D0*y1),maxm2(1))
!     call maxmod(maxm1(1),maxm2(1),mm1)
!   else if (mmopt==3) then ! Double minmod
!       call DOUBLEMINMOD(x1,y1,z1,mm1)
!   else 
!     call mc(2.0D0*x1,2.0D0*y1,z1,mm1)
!   end if
  
  HzjR(i,j)=hzt(i+2,j+1)+0.5D0*deta*mm1
  
  IF (j==1.AND.(CB(3)==4.OR.CB(3)==5)) THEN
    HzjR(i,1)=qA3(1,i)+zA3(i)
  END IF
  IF (j==Ny+1.AND.(CB(4)==4.OR.CB(4)==5)) THEN
    HzjL(i,Ny+1)=qA4(1,i)+zA4(i)
  END IF
  
End Do; End Do
!$omp end do

!3.Bathymethry Reconstruction using the linear reconstruction of the free surface

!---------loop 6--------------- no privates
!$omp do
Do i=1,Nx+1; Do j=1,Ny

  ZiL(i,j)=HziL(i,j)-QiL(1,i,j) !De la celda 1 a la N+1
  ZiR(i,j)=HziR(i,j)-QiR(1,i,j) !De la celda 0 a la N
  
  
  IF (i==1.AND.(CB(1)==4.OR.CB(1)==5)) THEN
  ZiR(1,j)=zA1(j)
  END IF
  IF (i==Nx+1.AND.(CB(2)==4.OR.CB(2)==5)) THEN
  ZiL(Nx+1,j)=zA2(j)
  END IF
End Do; End Do
!$omp end do NOWAIT

!---------loop 7---------------no privates
!$omp do
Do i=1,Nx; Do j=1,Ny+1
  ZjL(i,j)=HzjL(i,j)-QjL(1,i,j)
  ZjR(i,j)=HzjR(i,j)-QjR(1,i,j)
  
  IF (j==1.AND.(CB(3)==4.OR.CB(3)==5)) THEN
  ZjR(i,1)=zA3(i)
  END IF
  IF (j==Ny+1.AND.(CB(3)==4.OR.CB(4)==5)) THEN
  ZjL(i,Ny+1)=zA4(i)
  END IF
End Do; End Do
!$omp end do

!4.Hydrostatic Reconstruction of the water heigh

!hmas contiene h_i+1/2mas desde i=0+1/2 a i=N+1/2
!hmenos contiene h_i+1/2menos desde i=0+1/2 a i=N+1/2

!---------loop 8---------------zm
!$omp do
DO i=1,Nx+1; DO j=1,Ny
zm=max(ZiR(i,j),ZiL(i,j))	!Max(Z0R,Z1L)	!Max(ZNR,ZN+1L) !MAX(ZiR,Zi+1L)

hmenosi(i,j)=max(0.0D0,(QiR(1,i,j)+ZiR(i,j)-zm))	!max(0,h0R+Z0R-zm)	!max(0,hNR+ZNR-zm)	!hi+1/2- desde i=0 a i=Nx

hmasi(i,j)=max(0.0D0,(QiL(1,i,j)+ZiL(i,j)-zm))	!max(0,h1L+Z1L-zm)	!max(0,hN+1L+Zn+1L-zm) !hi-1/2+ desde i=1 a i=Nx+1

if (hmenosi(i,j)==0.0D0) then
QiR(2,i,j)=0.0D0
QiR(3,i,j)=0.0D0
end if

if (hmasi(i,j)==0.0D0) then
QiL(2,i,j)=0.0D0
QiL(3,i,j)=0.0D0
end if

END DO; END DO
!$omp end do NOWAIT

!---------loop 9---------------
!$omp do
DO i=1,Nx; DO j=1,Ny+1
zmj=max(ZjR(i,j),ZjL(i,j))

hmenosj(i,j)=max(0.0D0,(QjR(1,i,j)+ZjR(i,j)-zmj))

hmasj(i,j)=max(0.0D0,(QjL(1,i,j)+ZjL(i,j)-zmj))

if (hmenosj(i,j)==0.0D0) then
QjR(2,i,j)=0.0D0
QjR(3,i,j)=0.0D0
end if

if (hmasj(i,j)==0.0D0) then
QjL(2,i,j)=0.0D0
QjL(3,i,j)=0.0D0
end if
END DO; END DO
!$omp end do NOWAIT


!5.Fluxes Calculations: -Riemann solver using the values from the hydrostatic reconstruction
!			-Plus balancing source term


! !Metrics desde Xi0 a XiN+1 y Jacobianos
!Jac_global_xi: Jac desde Xi0+1/2 a Xi Nx+1/2 para j de 1 a Ny, Jacobiano en los bordes de la celda
!Jac_global_eta: Jac desde eta 0+1/2 a eta Ny+1/2 para i de 1 a Nx, Jacobiano en los bordes de la celda
!Jac_global: Jac de (xi,eta)=(1,1) a (Nx,Ny), Jacobiano centrado en la celda
! 

!---------loop 10---------------
!$omp do
DO i=2,Nx+2; DO j=3,Ny+2
xi(:,i-1,j-2)=xit(:,i,j)
!Jac1(i-1,j-2)=xit(1,i,j)*etat(2,i,j)-xit(2,i,j)*etat(1,i,j) !Ji+1/2,j, calculo J i=0 a i=Nx+1, j=1 a j=Ny
END DO; END DO
!$omp end do NOWAIT

!---------loop 11---------------
!$omp do
DO i=3,Nx+2; DO j=2,Ny+2
eta(:,i-2,j-1)=etat(:,i,j)
!Jac2(i-2,j-1)=xit(1,i,j)*etat(2,i,j)-xit(2,i,j)*etat(1,i,j) !Ji,j+1/2, calculo J j=0 a j=Ny+1, i=1 a i=Nx
END DO; END DO
!$omp end do


!Qstar
!Variables en cada Qstar, de la interfaz 0 a N

!---------loop 12--------------- qmenos,ximenos, qmas,ximas, qs
!$omp do
do i=1,Nx+1; do j=1,Ny
!Parto con la interfaz  0+1/2,j hasta la interfaz N+1/2, en total son N+1 interfaces

qmenos(1)=hmenosi(i,j)		!hmenosi contiene info desde la interfaz 0+1/2 hasta la N+1/2
qmenos(2)=QiR(2,i,j)		!QiR contiene info desde la celda 0 a N
qmenos(3)=QiR(3,i,j)
ximenos(:)=xi(:,i,j)		!xi contiene info desde la celda 0 a la N+1, aquí se usa de la 0 a la N (N+1)
qmas(1)=hmasi(i,j)		!hmasi contiene info desde la interfaz 0+1/2 hasta la N+1/2
qmas(2)=QiL(2,i,j)		!QiL contiene info desde la celda 1 a N+1
qmas(3)=QiL(3,i,j)
ximas(:)=xi(:,i+1,j)		!Aqui se usa la de la celda que sigue, de 1 a N+1 (N+1)

call vfroencv(i,j,FR2,qmas,qmenos,ximas,ximenos,hmin,qs)


Qstari(:,i,j)=qs(:)

! if (i==Nx+1.and.CB(2)==5) then
! Qstari(:,i,j)=qA2(:,j)
! end if

end do; end do
!$omp end do NOWAIT

!---------loop 13---------------qmenos,etamenos, qmas,etamas, qs
!$omp do
do i=1,Nx; do j=1,Ny+1

qmenos(1)=hmenosj(i,j)
qmenos(2)=QjR(2,i,j)
qmenos(3)=QjR(3,i,j)
etamenos(:)=eta(:,i,j)
qmas(1)=hmasj(i,j)
qmas(2)=QjL(2,i,j)
qmas(3)=QjL(3,i,j)
etamas(:)=eta(:,i,j+1)

call vfroencv(i,j,FR2,qmas,qmenos,etamas,etamenos,hmin,qs)

Qstarj(:,i,j)=qs(:)

end do; end do
!$omp end do
!FLUXES

!F fluxes

!---------loop 14---------------U1,Fs,Fb
!$omp do
DO i=1,Nx; DO j=1,Ny

!Fmas=F_i-1/2 + =F0+1!2+ a FN+1/2+

U1=Qstari(2,i,j)*xi(1,i,j)+Qstari(3,i,j)*xi(2,i,j)

Fs(1)=Qstari(1,i,j)*U1*1.0D0/Jac_global_xi(i,j)
Fs(2)=Qstari(1,i,j)*Qstari(2,i,j)*U1*1.0D0/Jac_global_xi(i,j)+&
0.5D0*(1.0D0/FR2)*(Qstari(1,i,j)**2.0D0)*xi(1,i,j)*1.0D0/Jac_global_xi(i,j)
Fs(3)=Qstari(1,i,j)*Qstari(3,i,j)*U1*1.0D0/Jac_global_xi(i,j)+&
0.5D0*(1.0D0/FR2)*(Qstari(1,i,j)**2.0D0)*xi(2,i,j)*1.0D0/Jac_global_xi(i,j)

Fb(1)=0.0D0
Fb(2)=1.0D0/(2.0D0*FR2)*(QiL(1,i,j)**2.0D0*xi(1,i+1,j)/&
Jac_global(i,j)-hmasi(i,j)**2.0D0*xi(1,i,j)*1.0D0/Jac_global_xi(i,j))
Fb(3)=1.0D0/(2.0D0*FR2)*(QiL(1,i,j)**2.0D0*xi(2,i+1,j)/&
Jac_global(i,j)-hmasi(i,j)**2.0D0*xi(2,i,j)*1.0D0/Jac_global_xi(i,j))


Fmas(:,i,j)=Fs+Fb
! print*,'Fstar',Qstari(1,i,j),U1
! print*,'holo', Fs, Fb, 'Jac',Jac_global(i,j), Jac_global_xi(i,j)
! pause

!Fmenos=F_i+1/2 -


U1=Qstari(2,i+1,j)*xi(1,i+1,j)+Qstari(3,i+1,j)*xi(2,i+1,j)

Fs(1)=Qstari(1,i+1,j)*U1*1.0D0/Jac_global_xi(i+1,j)
Fs(2)=Qstari(1,i+1,j)*Qstari(2,i+1,j)*U1*1.0D0/Jac_global_xi(i+1,j)+&
0.5D0*(1.0D0/FR2)*(Qstari(1,i+1,j)**2.0D0)*xi(1,i+1,j)*1.0D0/Jac_global_xi(i+1,j)
Fs(3)=Qstari(1,i+1,j)*Qstari(3,i+1,j)*U1*1.0D0/Jac_global_xi(i+1,j)+&
0.5D0*(1.0D0/FR2)*(Qstari(1,i+1,j)**2.0D0)*xi(2,i+1,j)*1.0D0/Jac_global_xi(i+1,j)

Fb(1)=0.0D0
Fb(2)=1.0D0/(2.0D0*FR2)*(QiR(1,i+1,j)**2.0D0*xi(1,i+1,j)/&
Jac_global(i,j)-hmenosi(i+1,j)**2.0D0*xi(1,i+1,j)/Jac_global_xi(i+1,j))
Fb(3)=1.0D0/(2.0D0*FR2)*(QiR(1,i+1,j)**2.0D0*xi(2,i+1,j)/&
Jac_global(i,j)-hmenosi(i+1,j)**2.0D0*xi(2,i+1,j)/Jac_global_xi(i+1,j))


Fmenos(:,i,j)=Fs+Fb

End do; End do
!$omp end do NOWAIT

!G fluxes

!---------loop 15---------------U2,Gs,Gb
!$omp do
DO i=1,Nx; DO j=1,Ny

U2=Qstarj(2,i,j)*eta(1,i,j)+Qstarj(3,i,j)*eta(2,i,j)

Gs(1)=Qstarj(1,i,j)*U2*1.0D0/Jac_global_eta(i,j)
Gs(2)=Qstarj(1,i,j)*Qstarj(2,i,j)*U2*1.0D0/Jac_global_eta(i,j)+&
0.5D0*(1.0D0/FR2)*(Qstarj(1,i,j)**2.0D0)*eta(1,i,j)*1.0D0/Jac_global_eta(i,j)
Gs(3)=Qstarj(1,i,j)*Qstarj(3,i,j)*U2*1.0D0/Jac_global_eta(i,j)+&
0.5D0*(1.0D0/FR2)*(Qstarj(1,i,j)**2.0D0)*eta(2,i,j)*1.0D0/Jac_global_eta(i,j)

Gb(1)=0.0D0
Gb(2)=1.0D0/(2.0D0*FR2)*(QjL(1,i,j)**2.0D0*eta(1,i,j+1)/&
Jac_global(i,j)-hmasj(i,j)**2.0D0*eta(1,i,j)/Jac_global_eta(i,j))
Gb(3)=1.0D0/(2.0D0*FR2)*(QjL(1,i,j)**2.0D0*eta(2,i,j+1)/&
Jac_global(i,j)-hmasj(i,j)**2.0D0*eta(2,i,j)/Jac_global_eta(i,j))

Gmas(:,i,j)=Gs+Gb

!Gmenos=G_j+1/2 --

U2=Qstarj(2,i,j+1)*eta(1,i,j+1)+Qstarj(3,i,j+1)*eta(2,i,j+1)


Gs(1)=Qstarj(1,i,j+1)*U2*1.0D0/Jac_global_eta(i,j+1)
Gs(2)=Qstarj(1,i,j+1)*Qstarj(2,i,j+1)*U2*1.0D0/Jac_global_eta(i,j+1)+&
  1.0D0/(2.0D0*FR2)*(Qstarj(1,i,j+1)**2.0D0)*eta(1,i,j+1)*1.0D0/Jac_global_eta(i,j+1)
Gs(3)=Qstarj(1,i,j+1)*Qstarj(3,i,j+1)*U2*1.0D0/Jac_global_eta(i,j+1)+&
  1.0D0/(2.0D0*FR2)*(Qstarj(1,i,j+1)**2.0D0)*eta(2,i,j+1)*1.0D0/Jac_global_eta(i,j+1)

Gb(1)=0.0D0
Gb(2)=1.0D0/(2.0D0*FR2)*(QjR(1,i,j+1)**2.0D0*eta(1,i,j+1)/Jac_global(i,j)-hmenosj(i,j+1)**2.0D0*eta(1,i,j+1)/Jac_global_eta(i,j+1))
Gb(3)=1.0D0/(2.0D0*FR2)*(QjR(1,i,j+1)**2.0D0*eta(2,i,j+1)/Jac_global(i,j)-hmenosj(i,j+1)**2.0D0*eta(2,i,j+1)/Jac_global_eta(i,j+1))

Gmenos(:,i,j)=Gs+Gb

END DO; END DO
!$omp end do NOWAIT

!6. Centered Source Term

!---------loop 16---------------no privates
!$omp do
DO i=1,Nx; DO j=1,Ny

!Balancing term for F fluxes
SF(1,i,j)=0.0D0
SF(2,i,j)=(-1.0D0/FR2)*xi(1,i+1,j)*(QiR(1,i+1,j)+QiL(1,i,j))*0.5D0*(ZiR(i+1,j)-ZiL(i,j))*(1.0D0/dxi)
SF(3,i,j)=(-1.0D0/FR2)*xi(2,i+1,j)*(QiR(1,i+1,j)+QiL(1,i,j))*0.5D0*(ZiR(i+1,j)-ZiL(i,j))*(1.0D0/dxi)
!Balancing term for G fluxes
SG(1,i,j)=0.0D0
SG(2,i,j)=(-1.0D0/FR2)*eta(1,i,j+1)*(QjR(1,i,j+1)+QjL(1,i,j))*0.5D0*(ZjR(i,j+1)-ZjL(i,j))*(1.0D0/deta)
SG(3,i,j)=(-1.0D0/FR2)*eta(2,i,j+1)*(QjR(1,i,j+1)+QjL(1,i,j))*0.5D0*(ZjR(i,j+1)-ZjL(i,j))*(1.0D0/deta)

!Full Balancing Source Term
SC(1,i,j)=0.0D0
SC(2,i,j)=SF(2,i,j)+SG(2,i,j)
SC(3,i,j)=SF(3,i,j)+SG(3,i,j)

END DO; END DO
!$omp end do

!$omp end parallel

END SUBROUTINE FLUXES


SUBROUTINE MINMOD(x,y,r)
real (kind=8) :: x, y, r
IF ((x>0.0D0).AND.(y>0.0D0)) THEN
  r=min(x,y)
ELSE IF ((x<0.0D0).AND.(y<0.0D0)) THEN
  r=max(x,y)
ELSE
  r=0.0D0
END IF
 
END SUBROUTINE MINMOD

SUBROUTINE MAXMOD(x,y,r)
real (kind=8) :: x, y, r
IF ((x>0.0D0).AND.(y>0.0D0)) THEN
  r=max(x,y)
ELSE IF ((x<0.0D0).AND.(y<0.0D0)) THEN
  r=min(x,y)
ELSE
  r=0.0D0
END IF
 
END SUBROUTINE MAXMOD

SUBROUTINE mc(x,y,z,r)
real (kind=8) :: x, y,z, r
IF ((x>0.0D0).AND.(y>0.0D0).AND.(z>0.0D0)) THEN
  r=min(x,y)
ELSE IF ((x<0.0D0).AND.(y<0.0D0).AND.(z<0.0D0)) THEN
  r=max(x,y)
ELSE
  r=0.0D0
END IF
 
END SUBROUTINE mc

SUBROUTINE DOUBLEMINMOD(x,y,z,r)
  real (kind=8) :: x, y, z, r, tmp1
  tmp1=min(2.0D0*x,2.0D0*y,z)
  r=max(0.0D0,tmp1)

END SUBROUTINE DOUBLEMINMOD

!1p 27.23
!2p 27.12