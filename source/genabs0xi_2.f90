! Funciones para el borde xi=0
!GA= 1:onda eta(t,x), 2: q(t,x)+h(t,x), 3: q(t,x)+u(t,x)
!Subrutina de la condicion de borde de generacion/absorcion de señales o hidrogramas
!2o y 4o paso RK
!Poner como argumento lo de fricci´on!

SUBROUTINE genabs0xi_1_2(fopt,tipo,MC,Nx,Ny,Fr2,dep,etaL,NL,h0,t,dt,q,zt,ep_x,qA0,qA,zA)
USE coords

!Señal, misma señal en todos los nodos, seña perpendicular al dominio y v=0.0, borde paralelo al eje y

!datos: eta(t), se asume igual para cada nodo del borde, 
!Aproximacion aguas someras
!etaL=señal
!q=qold sin celdas ficticias
!zt=cotas fondo con celdas ficticias
!ep_x=metricas sin celdas ficticias
!Borde xi=0

implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hL, uL, epL, epxL, RL,dep, zepL
real (kind=8)::epxA,zepA,tauL,epR,epxR,uR,hR
real (kind=8)::zepR, tauR,RR, Rmas, Rmenos,hAo,uAo,tauAo
integer:: Nx,Ny,j,i, NL, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA0,qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(NL,2)::etaL
real (kind=8), dimension(Nx,Ny)::MC

borde=1

! tL=epA+C0*dt
IF ((t+dt)>maxval(etaL(:,1))) THEN
    etai=0.0D0;
ELSE
    call interp1(NL,etaL(:,1),etaL(:,2),(t+dt),etai)
END IF
epA=0.0D0

DO j=1,Ny
DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO

!L
!h0=0.0D0-zt(3,j+2) !!!!SOLO MATAQUITO!!!!!
 C0=sqrt(h0/Fr2)
 
epL=epA-dt*C0
hL=etai+h0
uL=(hL/Fr2)**0.5D0*(etai)/hL
call interp1(Nx,coordxi,zep_x,epL,zepL)

!Situacion Borde
epx=ep_x(1,:,j)
call interp1(Nx,coordxi,epx,epA,epxA)
zepA=(zt(4,j+2)-zt(3,j+2))/dep

!RL
call interp1(Nx,coordxi,epx,epL,epxL)
RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL

!Rmas
! print*,'qA0=',qA0(1,j),qA0(2,j)
hAo=qA0(1,j)
uAo=qA0(2,j)
if (fopt==0) then
tauAo=0.0D0
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hAo,uAo,tauAo)
  tauAo=tauAo/hAo
  else
  tauAo=0.0D0
  end if

  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if


Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxL*tauL+epxA*tauAo)*dt

!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
!print*,h(1),j
!call fzero0_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)
call fzero_B2(borde,Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)

RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

if (fopt==0) then
tauR=0.0D0
else
  if (hR/=0.0D0) then
  C=MC(1,j)
  call  tau(tipo,C,Fr2,hR,uR,tauR)  
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if

Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxR*tauR+epxA*tauAo)*dt

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0

END DO


END SUBROUTINE genabs0xi_1_2
!----------------------------------------------------------
SUBROUTINE genabs0xi_2_2(fopt,tipo,MC,Nx,Ny,Fr2,dep,qs,hs,h0,Ns,t,dt,q,zt,ep_x,qA0,qA,zA)
!datos: q(x,t), h(x,t) para cada nodo del borde
!Borde eta=0
USE coords
implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, qL,hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,&
  tauL,epR,epxR,uR,hR,zepR, tauR,RR, Rmas, Rmenos,hAo,uAo,tauAo
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA0,qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(Ns,2)::qs,hs
real (kind=8), dimension(Nx,Ny)::MC
!allocate(q(3,Nx,Ny),qA(3,Ny),zt(Nx+4,Ny+4),zA(Ny),ep_x(2,Nx,Ny),h(Nx),u(Nx),etaL(NL,2),zep_x(Nx),epx(Nx))
borde=1

!Señal, misma señal en todos los nodos


! tL=epA+C0*dt
IF ((t+dt)>maxval(qs(:,1))) THEN
    print*, 'Faltan datos hidrograma'
    stop
ELSE
    call interp1(Ns,qs(:,1),qs(:,2),(t+dt),qL)
    call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hL)
END IF

 C0=sqrt(h0/Fr2)
epA=0.0D0
uL=qL/hL
epL=epA-dt*C0

DO j=1,Ny

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque zt incluye celdas ficticias
END DO

!L
call interp1(Nx,coordxi,zep_x,epL,zepL)


epx=ep_x(1,:,j)
!RL
call interp1(Nx,coordxi,epx,epL,epxL)
RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL

!Situacion Borde
call interp1(Nx,coordxi,epx,epA,epxA)

zepA=(zt(4,j+2)-zt(3,j+2))/dep

!Rmas
hAo=qA0(1,j)
uAo=qA0(2,j)
if (fopt==0) then
tauAo=0.0D0
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hAo,uAo,tauAo)
  tauAo=tauAo/hAo
  else
  tauAo=0.0D0
  end if

  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if

end if


Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxL*tauL+epxA*tauAo)*dt

!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
!!!call fzero0_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)
call fzero_B2(borde,Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)

RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

if (fopt==0) then
tauR=0.0D0
else
  if (hR/=0.0D0) then
  C=MC(1,j)
  call  tau(tipo,C,Fr2,hR,uR,tauR)  
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if

Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxL*tauR+epxA*tauAo)*dt

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0

END DO
END SUBROUTINE genabs0xi_2_2

!-------------------------------------------------------------------------------------------------------

!----------------------------------------------------------
SUBROUTINE genabs0xi_3_2(fopt,tipo,MC,Nx,Ny,Fr2,dep,us,hs,h0,Ns,t,dt,q,zt,ep_x,qA0,qA,zA)
!datos: q(x,t), h(x,t) para cada nodo del borde
!Borde eta=0
USE coords
implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, qL,hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,&
epR,epxR,uR,hR,zepR, tauR,RR, Rmas, Rmenos, hAo,uAo,tauAo
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA0,qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(Ns,2)::qs,hs, us
real (kind=8), dimension(Nx,Ny)::MC
!allocate(q(3,Nx,Ny),qA(3,Ny),zt(Nx+4,Ny+4),zA(Ny),ep_x(2,Nx,Ny),h(Nx),u(Nx),etaL(NL,2),zep_x(Nx),epx(Nx))
borde=1

!Señal, misma señal en todos los nodos


! tL=epA+C0*dt
IF ((t+dt)>maxval(us(:,1))) THEN
    print*, 'Faltan datos hidrograma'
    stop
!     uL=us(Ns,2)
!     hL=hs(Ns,2)
ELSE
    call interp1(Ns,us(:,1),us(:,2),(t+dt),uL)
    call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hL)
!     print*,'Ns=',Ns,uL,hL
END IF

 C0=sqrt(h0/Fr2)
epA=0.0D0

epL=epA-dt*C0

DO j=1,Ny

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque zt incluye celdas ficticias
END DO

!L
call interp1(Nx,coordxi,zep_x,epL,zepL)


epx=ep_x(1,:,j)
!RL
call interp1(Nx,coordxi,epx,epL,epxL)
RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL

!Situacion Borde
call interp1(Nx,coordxi,epx,epA,epxA)

zepA=(zt(4,j+2)-zt(3,j+2))/dep

!Rmas
hAo=qA0(1,j)
uAo=qA0(2,j)
if (fopt==0) then
tauAo=0.0D0
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hAo,uAo,tauAo)
  tauAo=tauAo/hAo
  else
  tauAo=0.0D0
  end if

  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if

end if


Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxL*tauL+epxA*tauAo)*dt

!Rmenos
h=q(1,:,j)	
u=q(2,:,j)

!print*, h(:)
!call fzero0_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)

call fzero_B2(Borde,Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)

RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

if (fopt==0) then
tauR=0.0D0
else
  if (hR/=0.0D0) then
  C=MC(1,j)
  call  tau(tipo,C,Fr2,hR,uR,tauR)  
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if

Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxL*tauR+epxA*tauAo)*dt

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0

END DO

! print*, qA(1,:)

END SUBROUTINE genabs0xi_3_2

! 
!-----------------------------------------------------------------------

SUBROUTINE genabs0xi_9_2(fopt,tipo,MC,Nx,Ny,Fr2,dep,etaL,timeS,NL,h0,t,dt,q,zt,ep_x,qA0,qA,zA)
USE coords
USE geometries
!Señal, misma señal en todos los nodos, seña perpendicular al dominio y v=0.0, borde paralelo al eje y

!datos: eta(t), se asume igual para cada nodo del borde, 
!Aproximacion aguas someras
!etaL=señal
!q=qold sin celdas ficticias
!zt=cotas fondo con celdas ficticias
!ep_x=metricas sin celdas ficticias
!Borde xi=0

implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,&
tauL,epR,epxR,uR,hR,zepR, tauR,RR, Rmas, Rmenos, hAo,uAo,tauAo
integer:: Nx,Ny,j,i, NL, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA0,qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(NL,2)::etaL, etaL_tmp
real (kind=8), dimension(NL,2)::timeS
real (kind=8), dimension(Nx,Ny)::MC

borde=1

DO j=1,Ny
! tL=epA+C0*dt
IF ((t+dt)>maxval(timeS)) THEN
    etai=0.0D0;
ELSE
    call interp1(NL,timeS,etaL(:,j),(t+dt),etai)
END IF
epA=0.0D0
!print*, h0
!pause
!DO j=1,Ny
DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO

!L
!h0=0.0D0-zt(3,j+2) !!!!SOLO MATAQUITO!!!!!
!ESSAI POUR LA CONDITION ENTRANTE!!!! CHANGEMENT DANS C0 et hL
!h0=h0-z_global(1,j)

 C0=sqrt((h0-z_global(1,j))/Fr2)
 
epL=epA-dt*C0
hL=etai+h0-z_global(1,j)
uL=(hL/Fr2)**0.5D0*(etai)/hL
call interp1(Nx,coordxi,zep_x,epL,zepL)

!Situacion Borde
epx=ep_x(1,:,j)
call interp1(Nx,coordxi,epx,epA,epxA)
zepA=(zt(4,j+2)-zt(3,j+2))/dep

!RL
call interp1(Nx,coordxi,epx,epL,epxL)
RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL

!Rmas
! print*,'qA0=',qA0(1,j),qA0(2,j)
hAo=qA0(1,j)
uAo=qA0(2,j)
if (fopt==0) then
tauAo=0.0D0
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hAo,uAo,tauAo)
  tauAo=tauAo/hAo
  else
  tauAo=0.0D0
  end if

  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if


Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxL*tauL+epxA*tauAo)*dt

!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
!print*,h(1),j
!!OLD FZERO call fzero0_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)
call fzero_B2(borde,Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)

RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

if (fopt==0) then
tauR=0.0D0
else
  if (hR/=0.0D0) then
  C=MC(1,j)
  call  tau(tipo,C,Fr2,hR,uR,tauR)  
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if

Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(epxR*tauR+epxA*tauAo)*dt

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0


END DO


END SUBROUTINE genabs0xi_9_2
!----------------------------------------------------------

!-------------------------------------------------------------------------------------------------------
! 
! SUBROUTINE genabs0xi_4_2(fopt,tipo,C,Nx,Ny,Fr2,dep,hs,h0,Ns,t,dt,q,zt,ep_x,qA0,qA,zA)
! USE coords
! 
! !ESTA CONDICION NO SE DEBE USAR, HAY QUE REVISARLA
! 
! !Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
! 
! implicit none
! real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hR, uR, epR, epxR, RR,dep, zepR, epxA,zepA,tauR,Rmenos, hb,hA, hAo,uAo,tauAo, zmin,uA, Rmas
! integer:: Nx,Ny,j,i, Ns, fopt, tipo
! real (kind=8), dimension(3,Nx,Ny)::q
! real (kind=8), dimension(Nx+4,Ny+4)::zt
! real (kind=8), dimension(2,Nx,Ny)::ep_x
! real (kind=8), dimension(3,Ny)::qA0,qA
! real (kind=8), dimension(Ny)::zA
! real (kind=8), dimension(Nx)::h,u,zep_x,epx
! real (kind=8), dimension(Ns,2)::hs
! 
! !Señal, misma señal en todos los nodos
! 
! IF ((t+dt)>maxval(hs(:,1))) THEN
! 
!      print*, 'Faltan datos hidrograma'
!      stop
! ELSE
!     call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hb)
! END IF
! 
! epA=0.0D0
! zmin=minval(zt(3,3:Ny+2))
! 
! DO j=1,Ny
! 
! epx=ep_x(1,:,j)
! call interp1(Nx,coordxi,epx,epA,epxA)
! 
! DO i=1,Nx
! zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
! END DO
! 
! !Rmenos
! h=q(1,:,j)	
! u=q(2,:,j)
! if (h(1) == 0.0D0) then  !Caso seco
! hR=0.0D0
! uR=0.0D0
! RR=0.0D0
! 
! hAo=qA0(1,j)
! uAo=qA0(2,j)
! else
!   !call fzero0_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)
!   call fzero_B2(Borde,Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epR,epxR,uR,hR,zepR)
!   RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR
! end if
! 
! if (fopt==0) then
! tauAo=0.0D0
! tauR=0.0D0
! else
!     if(hAo/=0.0D0) then
!     call tau(tipo,C,Fr2,hAo,uAo,tauAo)
!     tauAo=tauAo/hAo
!     else
!     tauAo=0.0D0
!     end if
!     
!     if(hR/=0.0D0) then
!     call  tau(tipo,C,Fr2,hR,uR,tauR)
!     tauR=tauR/hR
!     else
!     tauR=0.0D0
!     end if
!     
! end if
! 
! Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauR*epxR+tauAo*epxA)*dt
! 
! !Borde; 
! zA(j)=zt(3,j+2)
! 
! !Situacion Borde
! 
! if (hb.le.zA(j)) then
! hA=0.0D0
! uA=0.0D0
! else
! hA=(hb)-zA(j)
! Rmas=4.0D0*epxA*sqrt(hA/Fr2)+Rmenos
! uA=(Rmas+Rmenos)/(2.0D0*epxA)
! end if
! 
! qA(1,j)=hA
! qA(2,j)=uA
! qA(3,j)=0.0D0
! 
! END DO
! 
! 
! END SUBROUTINE genabs0xi_4_2
! 
! !-----------------------------------------------------------------------

