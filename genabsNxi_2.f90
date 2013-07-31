!MGP 20-04-11: FALTA REVISION!!! y CREAR LA GENABS9 y CAMBIAR LOS Fzero_B

! Funciones para el borde xi=N
!GA= 1:onda eta(t,x), 2: q(t,x)+h(t,x), 3: h(t,x)+u(t,x)
!Subrutina de la condicion de borde de generacion/absorcion de señales o hidrogramas
!2o y 4o paso RK


SUBROUTINE genabsNxi_1_2(fopt,tipo,MC,Nx,Ny,Fr2,dep,etaR,NR,h0,t,dt,q,zt,ep_x,qA0,qA,zA)
USE coords

!datos: eta(t), se asume igual para cada nodo del borde, 
!Aproximacion aguas someras
!etaL=señal
!q=qold sin celdas ficticias
!zt=cotas fondo con celdas ficticias
!ep_x=metricas sin celdas ficticias
!Borde xi=0

implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,epR,epxR,uR,hR,zepR, tauR,RR, Rmas, Rmenos, &
hAo,uAo,tauAo
integer:: Nx,Ny,j,i, NR, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA0,qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(NR,2)::etaR
real (kind=8), dimension(Nx,Ny)::MC
!Señal, misma señal en todos los nodos
borde=2

! tL=epA+C0*dt
IF ((t+dt)>maxval(etaR(:,1))) THEN
    etai=0.0D0;
ELSE
    call interp1(NR,etaR(:,1),etaR(:,2),(t+dt),etai)
END IF
epA=coordxi(Nx)+dep/2.0D0

DO j=1,Ny

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO
 C0=sqrt(h0/Fr2)

!R
hR=etai+h0
uR=(hR/Fr2)**0.5D0*(etai)/hR
epR=epA+dt*C0
call interp1(Nx,coordxi,zep_x,epR,zepR)

!Situacion Borde
epx=ep_x(1,:,j)
call interp1(Nx,coordxi,epx,epA,epxA)
zepA=(zt(Nx+2,j+2)-zt(Nx+1,j+2))/dep

!RR
call interp1(Nx,coordxi,epx,epR,epxR)
RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

!Rmenos
hAo=qA0(1,j)
uAo=qA0(2,j)
if (fopt==0) then
tauAo=0.0D0
tauR=0.0D0
else
    if(hAo/=0.0D0) then!Mojado
    C=MC(Nx,1)
    call tau(tipo,C,Fr2,hAo,uAo,tauAo)
    tauAo=tauAo/hAo
    else
    tauAo=0.0D0
    end if
    if (hR/=0.0D0) then
    C=MC(Nx,1)
    call tau(tipo,C,Fr2,hR,uR,tauR)
    tauR=tauR/hR
    else
    tauR=0.0D0
    end if
end if


Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauR*epxR+tauAo*epxA)*dt

!Rmas
h=q(1,:,j)	
u=q(2,:,j)
if (h(Nx) == 0.0D0) then  !Caso seco
hL=0.0D0
uL=0.0D0
RL=0.0D0
Rmas=0.0D0
else
!call fzeroN_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epL,epxL,uL,hL,zepL)
call fzero_B2(borde,Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epL,epxL,uL,hL,zepL)
RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL

if (fopt==0) then
tauL=0.0D0
else
    if(hL/=0.0D0) then
    C=MC(Nx,1)
    call  tau(tipo,C,Fr2,hL,uL,tauL)
    tauL=tauL/hL
    else
    tauL=0.0D0
    end if
end if

Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauL*epxL+tauAo*epxA)*dt

end if

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(Nx+3,j+2)+zt(Nx+2,j+2))/2.0D0

END DO


END SUBROUTINE genabsNxi_1_2
!----------------------------------------------------------
SUBROUTINE genabsNxi_2_2(fopt,tipo,MC,Nx,Ny,Fr2,dep,qs,hs,h0,Ns,t,dt,q,zt,ep_x,qA0,qA,zA)
USE coords
!USE Pich
!datos: qs(t), hs(t), se asume igual para cada nodo del borde, 
!Aproximacion aguas someras
!etaL=señal
!q=qold sin celdas ficticias
!zt=cotas fondo con celdas ficticias
!ep_x=metricas sin celdas ficticias
!Borde xi=0

implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,epR,epxR,qR,uR,hR,zepR, tauR,RR, Rmas, Rmenos, &
hAo,uAo,tauAo, zmin, AR, VR, Qt, dQ,hR1
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA0,qA
real (kind=8), dimension(Ny)::zA, haux
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(Ns,2)::qs,hs
real (kind=8), dimension(Nx,Ny)::MC

borde=2

!Señal, misma señal en todos los nodos


! tL=epA+C0*dt
IF ((t+dt)>maxval(qs(:,1))) THEN
    print*, 'Faltan datos hidrograma'
    stop
ELSE
    call interp1(Ns,qs(:,1),qs(:,2),(t+dt),qR)
    call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hR)
END IF

! !Pichilemu
! zmin=minval(zt(Nx+2,:))
! VR=qR/AR

  C0=sqrt(h0/Fr2)
epA=coordxi(Nx)+dep/2.0D0
epR=epA+dt*C0

DO j=1,Ny

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO

!R
call interp1(Nx,coordxi,zep_x,epR,zepR)
epx=ep_x(1,:,j)

!RR
call interp1(Nx,coordxi,epx,epR,epxR)
RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

!Situacion Borde
call interp1(Nx,coordxi,epx,epA,epxA)
zepA=(zt(Nx+2,j+2)-zt(Nx+1,j+2))/dep

!Rmenos
hAo=qA0(1,j)
uAo=qA0(2,j)

if (fopt==0) then
tauAo=0.0D0
tauR=0.0D0
else
    if(hAo/=0.0D0) then!Mojado
    C=MC(Nx,1)
    call tau(tipo,C,Fr2,hAo,uAo,tauAo)
    tauAo=tauAo/hAo
    else
    tauAo=0.0D0
    end if

    if(hR/=0.0D0) then!Mojado
    C=MC(Nx,1)
    call tau(tipo,C,Fr2,hR,uR,tauR)
    tauR=tauR/hR
    else
    tauR=0.0D0
    end if

end if

!Rmenos

if (hAo==0.0D0.AND.hR==0.0D0) then
Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)
!Rmenos=0.0D0
else if (hAo==0.0D0.AND.hR/=0.0D0) then
Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-(tauR*epxR)*dt

else if (hAo/=0.0D0.AND.hR==0.0D0) then
!Rmenos=0.0D0
Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauAo*epxA*dt
else
Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauR*epxR+tauAo*epxA)*dt
end if

!Rmas
h=q(1,:,j)	
u=q(2,:,j)
!print*, h(Nx),j

if (h(Nx) == 0.0D0) then  !Caso seco
hL=0.0D0
uL=0.0D0
RL=0.0D0
Rmas=0.0D0
else

  !call fzeroN_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epL,epxL,uL,hL,zepL)
  call fzero_B2(borde,Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epL,epxL,uL,hL,zepL)
  RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL
end if  

if (fopt==0) then
tauL=0
else
  if(hL/=0.0D0) then
  C=MC(Nx,1)
  call  tau(tipo,C,Fr2,hL,uL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if

!Rmas

if (hAo==0.0D0.AND.hL==0.0D0) then
Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)
!Rmas=0.0D0
else if (hAo==0.0D0.AND.hR/=0.0D0) then
Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-(tauL*epxL)*dt

else if (hAo/=0.0D0.AND.hR==0.0D0) then
Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-(tauAo*epxA)*dt
!Rmas=0.0D0
else
Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauL*epxL+tauAo*epxA)*dt

end if

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(Nx+3,j+2)+zt(Nx+2,j+2))/2.0D0

END DO

END SUBROUTINE genabsNxi_2_2

!----------------------------------------------------------
SUBROUTINE genabsNxi_3_2(fopt,tipo,MC,Nx,Ny,Fr2,dep,us,hs,h0,Ns,t,dt,q,zt,ep_x,qA0,qA,zA)
USE coords

!datos: us(t) hs(t), se asume igual para cada nodo del borde, 
!Aproximacion aguas someras
!etaL=señal
!q=qold sin celdas ficticias
!zt=cotas fondo con celdas ficticias
!ep_x=metricas sin celdas ficticias
!Borde xi=N

implicit none
real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA,qL, hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,epR,epxR,uR,hR,zepR, tauR,RR, Rmas, Rmenos,hAo,uAo,tauAo
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA, qA0
real (kind=8), dimension(Ny)::zA,haux
real (kind=8), dimension(Nx)::h,u,zep_x,epx
real (kind=8), dimension(Ns,2)::us,hs
real (kind=8), dimension(Nx,Ny)::MC
!allocate(q(3,Nx,Ny),qA(3,Ny),zt(Nx+4,Ny+4),zA(Ny),ep_x(2,Nx,Ny),h(Nx),u(Nx),etaL(NL,2),zep_x(Nx),epx(Nx))

borde=2
!Señal, misma señal en todos los nodos


IF ((t+dt)>maxval(us(:,1))) THEN
    print*, 'Faltan datos hidrograma'
    stop
ELSE
    call interp1(Ns,us(:,1),us(:,2),(t+dt),uR)
    call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hR)
END IF
 
 C0=sqrt(h0/Fr2)

epA=coordxi(Nx)+dep/2.0D0

epR=epA-dt*C0

DO j=1,Ny

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO

!R
call interp1(Nx,coordxi,zep_x,epR,zepR)

epx=ep_x(1,:,j)

!RR
call interp1(Nx,coordxi,epx,epR,epxR)
RR=uR*epxR-2.0D0*sqrt(hR/Fr2)*epxR

!Situacion Borde
call interp1(Nx,coordxi,epx,epA,epxA)
zepA=(zt(Nx+2,j+2)-zt(Nx+1,j+2))/dep

!Rmenos
hAo=qA0(1,j)
uAo=qA0(2,j)
if (fopt==0) then
tauAo=0.0D0
tauR=0.0D0
else
    if(hAo/=0.0D0) then!Mojado
    C=MC(Nx,1)
    call tau(tipo,C,Fr2,hAo,uAo,tauAo)
    tauAo=tauAo/hAo
    else
    tauAo=0.0D0
    end if
    if (hR/=0.0D0) then
    C=MC(Nx,1)
    call tau(tipo,C,Fr2,hR,uR,tauR)
    tauR=tauR/hR
    else
    tauR=0.0D0
    end if
end if

Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauR*epxR+tauAo*epxA)*dt

!Rmas
h=q(1,:,j)	
u=q(2,:,j)

if (h(Nx) == 0.0D0) then  !Caso seco
hL=0.0D0
uL=0.0D0
RL=0.0D0
Rmas=0.0D0
else
  
  !call fzeroN(Fr2,Nx,dep,epx,h,u,zep_x,dt,epL,epxL,uL,hL,zepL)
  call fzero_B2(borde,Fr2,Nx,dep,epx,h,u,zep_x,dt,epL,epxL,uL,hL,zepL)
  RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL
end if

if (fopt==0) then
tauL=0.0D0
else
    if(hL/=0.0D0) then
    C=MC(Nx,1)
    call  tau(tipo,C,Fr2,hL,uL,tauL)
    tauL=tauL/hL
    else
    tauL=0.0D0

    end if
end if

Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauL*epxL+tauAo*epxA)*dt

qA(1,j)=Fr2*(Rmas-Rmenos)**2.0D0/(16.0D0*epxA**2.0D0)
qA(2,j)=(Rmas+Rmenos)/(2.0D0*epxA)
qA(3,j)=0.0D0
zA(j)=(zt(Nx+3,j+2)+zt(Nx+2,j+2))/2.0D0

! Qt=Qt+qA(1,j)*qA(2,j)

END DO



END SUBROUTINE genabsNxi_3_2







! !-------------------------------------------------------------------------------------------------------
! SUBROUTINE genabsNxi_4_2(fopt,tipo,C,Nx,Ny,Fr2,dep,qs,h0,Ns,t,dt,q,zt,ep_x,qA0,qA,zA)
! USE coords
! !USE Pich
! !Compuerta aguas abajo, xi=N
! ! !ESTA CONDICION NO SE DEBE USAR, HAY QUE REVISARLA
! implicit none
! real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,Rmas, qb,hA, hAo,uAo,tauAo, zmin,uA
! integer:: Nx,Ny,j,i, Ns, fopt, tipo
! integer, dimension(1)::jmin
! real (kind=8), dimension(3,Nx,Ny)::q
! real (kind=8), dimension(Nx+4,Ny+4)::zt
! real (kind=8), dimension(2,Nx,Ny)::ep_x
! real (kind=8), dimension(3,Ny)::qA0,qA
! real (kind=8), dimension(Ny)::zA,Rcmas,hb
! real (kind=8), dimension(Nx)::h,u,zep_x,epx
! real (kind=8), dimension(Ns,2)::qs
! 
! !allocate(q(3,Nx,Ny),qA(3,Ny),zt(Nx+4,Ny+4),zA(Ny),ep_x(2,Nx,Ny),h(Nx),u(Nx),etaL(NL,2),zep_x(Nx),epx(Nx))
! 
! 
! !Señal, misma señal en todos los nodos
! 
! IF ((t+dt)>maxval(qs(:,1))) THEN
!     qb=qs(Ns,2)
! !     print*, 'Faltan datos hidrograma'
! !     stop
! ELSE
!     call interp1(Ns,qs(:,1),qs(:,2),(t+dt),qb)
! END IF
! !hR, Condición Hidráulica Compuerta Aguas Abajo
! 
! epA=coordxi(Nx)+dep/2.0D0
! zmin=minval(zt(Nx+2,3:Ny+2))
! jmin=minloc(zt(Nx+2,3:Ny+2))
! 
! DO j=1,Ny
! 
! epx=ep_x(1,:,j)
! call interp1(Nx,coordxi,epx,epA,epxA)
! 
! DO i=1,Nx
! zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
! END DO
! zepA=zep_x(Nx)
! 
! !Rmas
! h=q(1,:,j)	
! u=q(2,:,j)
! if (h(Nx) == 0.0D0) then  !Caso seco
! hL=0.0D0
! uL=0.0D0
! RL=0.0D0
! !Rmas=0.0D0
! hAo=qA0(1,j)
! uAo=qA0(2,j)
! else
!   call fzeroN_2(Fr2,Nx,dep,epx,h,u,zep_x,epxA,hAo,uAo,dt,epL,epxL,uL,hL,zepL)
!   !print*, uL,hL,epxL, zepL
!   RL=uL*epxL+2.0D0*sqrt(hL/Fr2)*epxL
! end if
! 
! if (fopt==0) then
! tauAo=0.0D0
! tauL=0.0D0
! else
!     if(hAo/=0.0D0) then
!     call tau(tipo,C,Fr2,hAo,uAo,tauAo)
!     tauAo=tauAo/hAo
!     else
!     tauAo=0.0D0
!     end if
!     
!     if(hL/=0.0D0) then
!     call  tau(tipo,C,Fr2,hL,uL,tauL)
!     tauL=tauL/hL
!     else
!     tauL=0.0D0
!     end if
!     
! end if
! 
! Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauL+tauAo)*dt*epxL
! 
! !Borde; Metricas y Zts y CI para ec. implicita
! Rcmas(j)=Rmas
! h0=q(1,Nx+2,jmin(1))
! zA(j)=zt(Nx+2,j+2)
! 
! !Situacion Borde
! if ((h0+zmin).le.zA(j)) then
! qA(1,j)=0.0D0
! qA(2,j)=0.0D0
! qA(3,j)=0.0D0
! else
! 
! call fzeroH(Rmas,qb,epxA,h0,zmin,zA(j),hA,uA)
! qA(1,j)=hA
! qA(2,j)=uA
! qA(3,j)=0.0D0
! end if
! 
! END DO
! 
! 
! END SUBROUTINE genabsNxi_4_2
! 
! !-------------------------------------------------------------------------------------------------------
! SUBROUTINE genabsNxi_5_2(fopt,tipo,C,Nx,Ny,Fr2,dep,qs,h0,Ns,t,dt,q,zt,ep_x,qA0,qA,zA)
! USE coords
! 
! !FUNCION solo para flujo unidimensional y xi_y=0.0
! !REVISAR!!!!
! 
! ! !ESTA CONDICION NO SE DEBE USAR, HAY QUE REVISARLA
! implicit none
! real (kind=8):: C,t,dt,Fr2,h0,C0,etai, epA, hL, uL, epL, epxL, RL,dep, zepL, epxA,zepA,tauL,Rmas, qb,hA, hAo,uAo,tauAo, zmin,uA, hminL,AL,TL,qL
! integer:: Nx,Ny,j,i, Ns, fopt, tipo
! integer, dimension(1)::jmin
! real (kind=8), dimension(3,Nx,Ny)::q
! real (kind=8), dimension(Nx+4,Ny+4)::zt
! real (kind=8), dimension(2,Nx,Ny)::ep_x
! real (kind=8), dimension(3,Ny)::qA0,qA
! real (kind=8), dimension(Ny)::zA,Rcmas,hb
! real (kind=8), dimension(Nx)::h,u,zep_x,epx,ep
! real (kind=8), dimension(Ns,2)::qs
! 
! !allocate(q(3,Nx,Ny),qA(3,Ny),zt(Nx+4,Ny+4),zA(Ny),ep_x(2,Nx,Ny),h(Nx),u(Nx),etaL(NL,2),zep_x(Nx),epx(Nx))
! 
! 
! !Señal, misma señal en todos los nodos
! 
! IF ((t+dt)>maxval(qs(:,1))) THEN
!     qb=qs(Ns,2)
! !     print*, 'Faltan datos hidrograma'
! !     stop
! ELSE
!     call interp1(Ns,qs(:,1),qs(:,2),(t+dt),qb)
! END IF
! !hR, Condición Hidráulica Compuerta Aguas Abajo
! 
! epA=coordxi(Nx)+dep/2.0D0
! zmin=minval(zt(Nx+2,3:Ny+2))
! jmin=minloc(zt(Nx+2,3:Ny+2))
! hminL=q(1,Nx+2,jmin(1))
! 
! call Saliente(zmin,jmin(1),hminL,hL,AL,TL,qL)
! 
! uL=qL/AL
! hL=AL/TL
! epxL=coordxi(Nx)
! RL=uL*epxL+2.0D0*sqrt(AL/TL/Fr2)*epxL
! 
! DO i=1,Nx
! zep_x(i)=(zt(i+3,jmin(1))-zt(i+2,jmin(1)))/dep	!j+2 porque incluye celdas ficticias
! END DO
! zepA=zep_x(Nx)
! 
! ep(1)=0.5D0
! DO i=2,Nx
! ep(i)=ep(i-1)+dep
! END DO
! epL=ep(Nx)-0.5D0
! call interp1(Nx,ep,zep_x,epL,zepL)
! 
! epx=ep_x(1,:,j)
! call interp1(Nx,coordxi,epx,epA,epxA)
! 
! 
! if (fopt==0) then
! tauAo=0.0D0
! tauL=0.0D0
! else
!     if(hAo/=0.0D0) then
!     call tau(tipo,C,Fr2,hAo,uAo,tauAo)
!     tauAo=tauAo/hAo
!     else
!     tauAo=0.0D0
!     end if
!     
!     if(hL/=0.0D0) then
!     call  tau(tipo,C,Fr2,hL,uL,tauL)
!     tauL=tauL/hL
!     else
!     tauL=0.0D0
!     end if
!     
! end if
! 
! Rmas=RL-0.5D0*dt/Fr2*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-0.5D0*(tauL+tauAo)*dt*epxL
! 
! h0=q(1,Nx+2,jmin(1))
! 
! if(qb/=0.0D0) then
! 
! call fzeroH_2(Rmas,qb,epxA,h0,zmin,hA,uA)
! 
! else
! hA=h0
! uA=uL
! end if
! 
! DO j=1,Ny
! zA(j)=zt(Nx+2,j+2)
! !Situacion Borde
! !Informacion entrante
! if (zA(j)<(hA+zmin)) then !Nodo mojado
! qA(1,j)=hA+zmin-zA(j)
! qA(2,j)=uA
! qA(3,j)=0.0D0
! else
! qA(1,j)=0.0D0
! qA(2,j)=0.0D0
! qA(3,j)=0.0D0
! end if
! END DO
! 
! 
! 
! END SUBROUTINE genabsNxi_5_2