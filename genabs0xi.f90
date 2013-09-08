! Funciones para el borde xi=0
!GA= 1:onda eta(t,x), 2: q(t,x)+h(t,x), 3: q(t,x)+u(t,x)
!Subrutina de la condicion de borde de generacion/absorcion de señales o hidrogramas
!1er y 3er paso RK
!Poner como argumento lo de fricci´on!

SUBROUTINE genabs0xi_1_1(fopt,tipo,MC,Nx,Ny,Fr2,dep,etaL,NL,h0,t,dt,q,zt,ep_x,qA,zA)
USE coords

!Señal, misma señal en todos los nodos, seña perpendicular al dominio y v=0.0, borde paralelo al eje y
! xi_y=0.0D0

!datos: eta(t), se asume igual para cada nodo del borde, 
!Aproximacion aguas someras
!etaL=señal
!q=qold sin celdas ficticias
!zt=cotas fondo con celdas ficticias
!ep_x=metricas sin celdas ficticias
!Borde xi=0

implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,epR,epxR,uR,hR,zepR, tauR,RR, Rmas, Rmenos
integer:: Nx,Ny,j,i, NL, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(NL,2)::etaL
real (kind=8), dimension(Nx,Ny)::MC
borde=1
!Señal, misma señal en todos los nodos
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

!h0=0.0D0-zt(3,j+2)!!!!SOLO MATAQUITO!!!!!

 C0=sqrt(h0/Fr2)
!L
epL=epA-dt*C0
hL=etai+h0
uL=sqrt(hL/Fr2)*(etai)/hL
call interp1(Nx,coordxi,zep_x,epL,zepL)

!Situacion Borde
epx=ep_x(1,:,j)
call interp1(Nx,coordxi,epx,epA,epxA)
zepA=(zt(4,j+2)-zt(3,j+2))/dep

!RL
call interp1(Nx,coordxi,epx,epL,epxL)
RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL


!Rmas
if (fopt==0) then
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if


Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*epxL*dt

!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
! print*, etai,h(1),j

!call fzero0(Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)
call fzero_B(borde,Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)

RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

if (fopt==0) then
tauR=0
else

  if (hR/=0.0D0) then
  C=MC(1,j)
  call  tau(tipo,C,Fr2,hR,uR,tauR)  
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if

Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*dt*epxR

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0

zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0

END DO


END SUBROUTINE genabs0xi_1_1
!----------------------------------------------------------
SUBROUTINE genabs0xi_2_1(fopt,tipo,MC,Nx,Ny,Fr2,dep,qs,hs,h0,Ns,t,dt,q,zt,ep_x,qA,zA)
!datos: q(x,t), h(x,t) para cada nodo del borde
!Borde eta=0
USE coords
implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA,qL,hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,epR,epxR,uR,hR,zepR,tauR, RR, Rmas, Rmenos
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(Ns,2)::qs,hs
real (kind=8), dimension(Nx,Ny)::MC
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
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
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
if (fopt==0) then
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if

Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*epxL*dt

!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
!call fzero0(Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)
call fzero_B(borde,Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)

RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR
if (fopt==0) then
tauR=0
else
  if (hR/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hR,uR,tauR)
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if

Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*epxR*dt

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0

END DO
 
END SUBROUTINE genabs0xi_2_1

!----------------------------------------------------------
!----------------------------------------------------------
SUBROUTINE genabs0xi_3_1(fopt,tipo,MC,Nx,Ny,Fr2,dep,us,hs,h0,Ns,t,dt,q,zt,ep_x,qA,zA)
!datos: u(y,t), h(y,t) para cada nodo del borde
!Borde eta=0
USE coords
implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA,qL,hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,epR,epxR,uR,hR,zepR,tauR, RR, Rmas, Rmenos
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(Ns,2)::qs,hs, us
real (kind=8), dimension(Nx,Ny)::MC
!Señal, misma señal en todos los nodos

borde=1
! tL=epA+C0*dt
IF ((t+dt)>maxval(us(:,1))) THEN
    print*, 'Faltan datos hidrograma'
    stop
!     uL=us(Ns,2)
!     hL=hs(Ns,2)
    
ELSE
    call interp1(Ns,us(:,1),us(:,2),(t+dt),uL)
    call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hL)
!     print*, uL,hL
!     pause
END IF
 
 
 
 C0=sqrt(h0/Fr2)
epA=0.0D0

epL=epA-dt*C0

DO j=1,Ny

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
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
if (fopt==0) then
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if

Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*epxL*dt
! print*, 'holi',j
! print*, hL,uL,tauL

!Rmenos
h=q(1,:,j)	
u=q(2,:,j)

!call fzero0(Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)
call fzero_B(borde,Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)

!print*, uR, hR,epR
RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

if (fopt==0) then
tauR=0.0D0
else
  if (hR/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hR,uR,tauR)
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if

Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*epxR*dt

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0
!print*, Rmas, Rmenos, epxA
!print*, qA(1,j), qA(2,j), j
END DO


! print*, qA(1,:), qA(2,:)
! pause

END SUBROUTINE genabs0xi_3_1


!------------------------------------------------------------
! Leandro Suarez BC
SUBROUTINE genabs0xi_9_1(fopt,tipo,MC,Nx,Ny,Fr2,dep,etaL,timeS,NL,h0,t,dt,q,zt,ep_x,qA,zA)
USE coords
USE geometries

!Señal, distinta señal en cada nodo, seña perpendicular al dominio y v=0.0, borde paralelo al eje y
! xi_y=0.0D0

!datos: eta(t), diferente para cada punto de entrada
!Aproximacion aguas someras
!etaL=señal
!q=qold sin celdas ficticias
!zt=cotas fondo con celdas ficticias
!ep_x=metricas sin celdas ficticias
!Borde xi=0

implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai,epA,hL,uL,epL,epxL,RL,dep,zepL,epxA,zepA,tauL,epR,epxR,uR,hR,zepR,tauR,RR,Rmas,Rmenos
integer:: Nx,Ny,j,i, NL, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(NL,Ny)::etaL
real (kind=8), dimension(NL)::timeS
real (kind=8), dimension(Nx,Ny)::MC
borde=1
!Señal, misma señal en todos los nodos
! tL=epA+C0*dt
!print*,'etaL', etaL(:,1)
!print*,'etaL2', etaL(:,2)

DO j=1,Ny

  IF ((t+dt)>maxval(timeS)) THEN
      etai=0.0D0;    
  ELSE
  ! 	print*,'NL', NL
  ! 	print*,'t+dt', t+dt
  ! 	print*,'etai', etai
  ! 	print*,'etaL',shape(etaL(:,j))
  ! 	pause
      call interp1(NL,timeS,etaL(:,j),(t+dt),etai)
      
      !print*,'etai', etai
  END IF
  epA=0.0D0
!   DO j=1,Ny
  DO i=1,Nx
  zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
  END DO

  !ESSAI POUR LA CONDITION ENTRANTE!!!!
  !h0=h0-z_global(1,j)

  C0=sqrt((h0-z_global(1,j))/Fr2)

  epL=epA-dt*C0
  hL=etai+h0-z_global(1,j)
  !si h0=0 y -zglobal(1,j) es la profundidad del fondo estamos bien
  !el problema es cuando eta no esta 'centrado' en el 0 del reposo del mar...
  !ahi habria que editar READGA
  !h0=sealevel
  uL=sqrt(hL/Fr2)*(etai)/hL
  call interp1(Nx,coordxi,zep_x,epL,zepL)

  !Situacion Borde
  epx=ep_x(1,:,j)
  call interp1(Nx,coordxi,epx,epA,epxA)
  zepA=(zt(4,j+2)-zt(3,j+2))/dep

  !RL
  call interp1(Nx,coordxi,epx,epL,epxL)
  RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL


  !Rmas
  if (fopt==0) then
  tauL=0.0D0
  else
    if (hL/=0.0D0) then
    C=MC(1,j)
    call tau(tipo,C,Fr2,hL,uL,tauL)
    tauL=tauL/hL
    else
    tauL=0.0D0
    end if
  end if



  Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*epxL*dt

  !Rmenos
  h=q(1,:,j)	
  u=q(2,:,j)
  ! print*, etai,h(1),j

  !call fzero0(Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)
  call fzero_B(borde,Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR) ! essai méthode de Brent
  RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

  if (fopt==0) then
    tauR=0
  else

    if (hR/=0.0D0) then
    C=MC(1,j)
    call  tau(tipo,C,Fr2,hR,uR,tauR)  
    tauR=tauR/hR
    else
    tauR=0.0D0
    end if
  end if


  Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*dt*epxR

  qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
  qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
  qA(3,j)=0.0D0
  !print*,(zt(3,j+2)+zt(2,j+2))/2.0D0
  zA(j)=(zt(3,j+2)+zt(2,j+2))/2.0D0
  !print*,j
END DO


END SUBROUTINE genabs0xi_9_1
!----------------------------------------------------------

! !----------------------------------------------------------
! 
! SUBROUTINE genabs0xi_4_1(fopt,tipo,C,Nx,Ny,Fr2,dep,hs,h0,Ns,t,dt,q,zt,ep_x,qA,zA)
! USE coords
! !ESTA CONDICION NO SE DEBE USAR, HAY QUE REVISARLA
! !Nivel fijo aguas afuera para flujo unidireccional v=0.0D0
! 
! implicit none
! real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hR, uR, epR, epxR, RR,dep, zepR, epxA,zepA,tauR,Rmenos, hA,zmin, uA, hb, Rmas
! integer:: Nx,Ny,j,i, Ns, fopt, tipo
! real (kind=8), dimension(3,Nx,Ny)::q
! real (kind=8), dimension(Nx+4,Ny+4)::zt
! real (kind=8), dimension(2,Nx,Ny)::ep_x
! real (kind=8), dimension(3,Ny)::qA
! real (kind=8), dimension(Ny)::zA, Rcmenos
! real (kind=8), dimension(Nx)::h,u,zep_x,epx
! real (kind=8), dimension(Ns,2)::hs
! 
! !Busco Qb
! IF ((t+dt)>maxval(hs(:,1))) THEN
!      print*, 'Faltan datos hidrograma'
!      stop
! ELSE
!     call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hb)
!     !hb es la cota (hs+zmin)
!     END IF
! 
! 
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
! ! Caracteristica Saliente
! !Rmenos
! h=q(1,:,j)	
! u=q(2,:,j)
! if (h(1) == 0.0D0) then  !Caso seco
! hR=0.0D0
! uR=0.0D0
! RR=0.0D0
! !Rmas=0.0D0
! else
!   call fzero0(Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)
!   RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR
! end if
!      
! if (fopt==0) then
! tauR=0.0D0
! else
!     if(hR/=0.0D0) then
!     call  tau(tipo,C,Fr2,hR,uR,tauR)
!     Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR/hR*dt*epxR
!     else
!     tauR=0.0D0
!     Rmenos=RR-0.5D0*dt/Fr2*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)
!     end if
! end if
! 
! 
! !Borde: Metricas y Zts y CI para ec. implicita
! Rcmenos(j)=Rmenos
! h0=h(1)
! zA(j)=zt(3,j+2)
! 
! !Situacion Borde
! !Informacion entrante: altura fija.
! 
! if ((hb).le.zA(j)) then !nodo seco
! hA=0.0D0
! uA=0.0D0
! else
! hA=(hb)-zA(j)
! Rmas=4.0D0*epxA*sqrt(hA/Fr2)+Rmenos
! uA=(Rmas+Rmenos)/(2.0D0*epxA)
! end if
! qA(1,j)=hA
! qA(2,j)=uA
! qA(3,j)=0.0D0
! 
! 
! END DO
! 
! 
! 
! END SUBROUTINE genabs0xi_4_1
! 
