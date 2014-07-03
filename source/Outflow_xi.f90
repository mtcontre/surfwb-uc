SUBROUTINE OUTFLOW0_xi(pasoRK,fopt,tipo,MC,Nx,Ny,Fr2,dep,etas1,timeS,Ns,t,dt,q,zt,ep_x,ep2_x,qA,zA)
!Basado en la rutina de Brett Sanders

!Subrutina que ingresa un q(t) unitario en cada nodo, en el borde xi=0, 
!Fijarse en el signo del caudal!!!
!Si q(t) entra qx o qy deben tener signo mas
!Si q(t) sale qx o qy deben apuntar hacia afuera, signo neg

USE coords
USE time0
!Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
implicit none
real (kind=8):: C,t,dt,Fr2,epA, hR, uR,vR, epR, epxR,epyR, RR,dep, &
zepR, epxA,epyA,epxA2,epyA2,zepA,tauR,Rmenos,hA,zmin,uA,vA, &
Rmas,RmenosU,U1, alfa,a1,a2,a3,a4,K,qp,U1u,uU,vU,hU, cero, JacA, us, vs, hs
real (kind=8)::etaL,uL,vL,hL,RL, epL,epxL,epyL, zepL,C0, h0, tauL, qpar, Uper,Upar, UparL, UparR 

integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde, pasoRK
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x, ep2_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA,qx,qy, eta1
real (kind=8), dimension(Nx)::h,u,v,zep_x,epx,epy, epx2, epy2
real (kind=8), dimension(Ny,Ns)::etas1
real (kind=8), dimension(Ns)::timeS,etas
real (kind=8), dimension(Nx,Ny)::MC
real (kind=8):: zumbral, x1,x2, x3, qx2, qy2, u0,v0

borde=1
zumbral=0.0D0+39.309606D0!38.848857D0!39.837508D0!39.899902D0 !206.5124D0!208.30362D0 !255.6058D0!326.0532D0
!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(timeS)) THEN
!      hb=hs(Ns,2)
!      print*, 'Faltan datos hidrograma'
!      stop
    Do j=1,Ny
    eta1(j)=1.5D0
    end do
    

ELSE

    !En el caso de que se tenga un qx y qy para cada nodo del borde
    !esto se debe hacer para cada j, usando
    Do j=1,Ny
    etas=etas1(j,:)
    call interp1(Ns,timeS(:),etas(:),(t+dt),eta1(j))
    end Do
END IF

!Borde

epA=0.0D0

zmin=minval(zt(3,3:Ny+2))

DO j=1,Ny

!Borde
zA(j)=zt(3,j+2)
h0=abs(zA(j)-39.309606D0)!+1.5D0!38.848857D0!39.881438D0)+1.5D0!-39.837508D0)+1.5D0!-206.5124D0)+1.4D0!208.30362D0)+1.4D0!) !Si es high tide
etaL=eta1(j)

epx=ep_x(1,:,j)
epy=ep_x(2,:,j)

epx2=ep2_x(1,:,j)
epy2=ep2_x(2,:,j)

call interp1(Nx,coordxi,epx,epA,epxA) !Xi_x
call interp1(Nx,coordxi,epy,epA,epyA) !Xi_y

call interp1(Nx,coordxi,epx2,epA,epxA2) !Eta_x
call interp1(Nx,coordxi,epy2,epA,epyA2) !Eta_y

JacA=epxA*epyA2-epyA*epxA2 

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO

zepA=(zt(3,j+2)-zt(2,j+2))/dep
!zepA=(zt(3,j+2)-zt(2,j+2))/dep 

!-----------------------------
! Caracteristica Saliente
!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
v=q(3,:,j)

!Solo tsunami
if (zA(j)>=zumbral) then !Borde Abierto
qA(1,j)=h(1)
qA(2,j)=u(1)
qA(3,j)=v(1)

else


if (h(1) == 0.0D0) then  !Caso seco
hR=0.0D0
uR=0.0D0
vR=0.0D0
RR=0.0D0

else

!   call fzeroEP(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epR,epxR,epyR,uR,vR,hR,zepR)
  !Cambiar fzeroEP a la perpendicular al borde y resolver en el plano real
  call fzeroEP(j,borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epR,epxR,epyR,uR,vR,hR,zepR)
  !call fzeroEP_brent(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epR,epxR,epyR,uR,vR,hR,zepR)
  
!   ! Usando la direccion Xi como la perpendicular al borde
    RR=(uR*epxR+vR*epyR)-2.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))
   
end if
     
if (fopt==0) then
tauR=0.0D0
Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
else
    if(hR/=0.0D0) then
    C=MC(1,j)
    call  tauU(tipo,C,Fr2,hR,uR,vR,tauR)
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+&
    zepA*(epxA**2.0D0+epyA**2.0D0))-tauR/hR*dt*(epxR+epyR)
    else
    tauR=0.0D0
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if

end if

! Caracteristica Entrante

 C0=sqrt(h0/Fr2)
epL=epA-dt*C0
hL=etaL+h0
!print*, j, h0, hL, h(1)

call interp1(Nx,coordxi,zep_x,epL,zepL)
call interp1(Nx,coordxi,epx,epL,epxL)
call interp1(Nx,coordxi,epy,epL,epyL)

! ---------------------------------------
! Usando metodología de Brett Sanders
! ---------------------------------------
!Metodologia Alturas
!-----------------------------------------
! Completo con hR, uR, vR el q0 en la primera iteracion
call q0(borde,pasoRK,epxA,epyA,hR,uR,vR,j)

!Solucion Ec. Orden 3
hU=q0_global(1,1,j)
uU=q0_global(2,1,j)
vU=q0_global(3,1,j)
 
!U perpendicular U1 asumiendo que Uper=U1 en la direccion Xi
U1u=uU*epxR0+vU*epyR0
RmenosU=U1u-2.0D0*sqrt(hU/Fr2)*sqrt(epxR0**2.0D0+epyR0**2.0D0)

RL=RmenosU+4.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))

if (fopt==0) then
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tau(tipo,C,Fr2,hL,uL,vL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if
Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*(epxL+epyL)*dt

!!!---------------------------------------------------------
! Solucion:
!En la direccion perpendicular dominio curvilíneo Xi
hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)
! En la direccion perpendicular dominio real
!hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0
Uper=(Rmas+Rmenos)/2.0D0
! print*,'Rmas=', Rmas, 'Rmenos=',Rmenos, 'hA=',hA,'U1A=',U1
! pause
! print*,'Uper=',U1

UparL=uL*epxA2+vL*epyA2 !U2
UparR=uR*epxA2+vR*epyA2 !U2
!UparL=-uL*sin(alfa)+vL*cos(alfa) !Esta debiese coincidir con U2
!UparR=-uR*sin(alfa)+vR*cos(alfa) !Esta debiese coincidir con U2

if (Uper>0.0D0) then
Upar=UparL
! print*,'Upar=',Upar
else
Upar=UparR
! print*,'Upar=',Upar
end if

uA=Uper*epyA2/JacA+Upar*-epyA/JacA !U1*eta_y/J-Upar*Xiy/Jac
vA=Uper*-epxA2/JacA+Upar*epxA/JacA

if (abs(uA)<10e-4) then
uA=0.0D0
end if
if (abs(vA)<10e-4) then
vA=0.0D0
end if



if (hA<1e-7) then
qA(1,j)=0.0D0
qA(2,j)=0.0D0
qA(3,j)=0.0D0
else
qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA
end if




!Esto es para la modelacion de tsunami solamente!!



end if
! print*,qA(:,j)
! pause
! ! ! ! !
! ! print*,j
! print*,'hL=',hL,'hR=',hR
! print*, 'hA=',hA, 'uA=',uA, 'vA=',vA
! !print*, 'qxA=',hA*uA, 'qyA=',hA*vA!, 'dh0=',hA-h0,'dhR=',hA-hR
! ! print*, 'qx=',qx(j), 'qy=',qy(j)!, 'etaL=',etaL
! ! 
!  print*, 'dh0=',hA-h0,'dhR=',hA-hR, 'dhL=',hL-hA
! print*, 'etaL=',etaL
! 
! pause

END DO




END SUBROUTINE OUTFLOW0_xi


! SUBROUTINE OUTFLOW0_xi(fopt,tipo,MC,Nx,Ny,Fr2,dep,hs,Ns,t,dt,q,zt,ep_x,qA,zA)
! !Subrutina que fija una altura de agua en el borde 0 del dominio
! 
! !Crear una que fije la altura en el borde N del dominio
! 
! USE coords
! USE time0
! !Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
! 
! implicit none
! real (kind=8):: C,t,dt,Fr2,epA, hR, uR,vR, epR, epxR,epyR, RR,dep, zepR, epxA,epyA,zepA,tauR,Rmenos, hb,hA,zmin,uA,vA, Rmas, U1, alfa,hU,uU,vU,U1u,RmenosU, uL,vL,UparL,UparR,Upar
! integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
! real (kind=8), dimension(3,Nx,Ny)::q
! real (kind=8), dimension(Nx+4,Ny+4)::zt
! real (kind=8), dimension(2,Nx,Ny)::ep_x
! real (kind=8), dimension(3,Ny)::qA
! real (kind=8), dimension(Ny)::zA
! real (kind=8), dimension(Nx)::h,u,v,zep_x,epx,epy
! real (kind=8), dimension(Ns,2)::hs
! real (kind=8), dimension(Nx,Ny)::MC
! borde=1
! !Busco Altura fija, se entrega la cota de la superficie libre
! 
! IF ((t+dt)>maxval(hs(:,1))) THEN
! !      hb=hs(Ns,2)
!      print*, 'Faltan datos hidrograma'
!      stop
! ELSE
!     call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hb)
!     !hb es la profundidad (hs)
!     !hs=(t,hb)
! END IF
! epA=0.0D0
! zmin=minval(zt(3,3:Ny+2))
! 
! DO j=1,Ny
! 
! 
! zA(j)=zt(3,j+2)
! 
! if ((hb).le.zA(j)) then !nodo seco
! hA=0.0D0
! uA=0.0D0
! vA=0.0D0
! 
! else !Calculo la info entrante 
! 
! epx=ep_x(1,:,j)
! epy=ep_x(2,:,j)
! call interp1(Nx,coordxi,epx,epA,epxA)
! call interp1(Nx,coordxi,epy,epA,epyA)
! 
! DO i=1,Nx
! zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
! END DO
! 
! ! Caracteristica Saliente
! !Rmenos
! h=q(1,:,j)	
! u=q(2,:,j)
! v=q(3,:,j)
! if (h(1) == 0.0D0) then  !Caso seco
! hR=0.0D0
! uR=0.0D0
! vR=0.0D0
! RR=0.0D0
! !Rmas=0.0D0
! else
!   !print*, j
!   !call fzeroEP(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epR,epxR,epyR,uR,vR,hR,zepR)
!   call fzeroEP_brent(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epR,epxR,epyR,uR,vR,hR,zepR)
!   RR=(uR*epxR+vR*epyR)-2.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))
! end if
!      
! if (fopt==0) then
! tauR=0.0D0
! Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
! else
!     if(hR/=0.0D0) then
!     C=MC(1,j)
!     call  tauU(tipo,C,Fr2,hR,uR,vR,tauR)
!     Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))-tauR/hR*dt*(epxR+epyR)
!     else
!     tauR=0.0D0
!     Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
!     end if
! end if
! 
! alfa=angulo1(j)
! 
! !Situacion Borde
! !Informacion entrante: altura fija.
! !Info en t=0
! hU=q0_global(1,1,j)
! uU=q0_global(2,1,j)
! vU=q0_global(3,1,j)
! 
! U1u=uU*epxA+vU*epyA
! 
! RmenosU=U1u-2.0D0*sqrt(hU/Fr2)*(epxA**2.0D0+epyA**2.0D0)
! 
! 
! !hA=(hb)-zA(j)
! Rmas=RmenosU+4.0D0*sqrt(epxA**2.0D0+epyA**2.0D0)*sqrt(hb/Fr2)
! hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)
! U1=(Rmas+Rmenos)/2.0D0
! 
! 
! uL=U1*cos(alfa)
! vL=U1*sin(alfa)
! 
! UparL=-uL*sin(alfa)+vL*cos(alfa)
! UparR=-uR*sin(alfa)+vR*cos(alfa)
! 
! if (U1>0.0D0) then
! Upar=UparL
! else
! Upar=UparR
! end if
! 
! uA=U1*cos(alfa)-Upar*sin(alfa)
! vA=U1*sin(alfa)+Upar*cos(alfa)
! 
! if (abs(uA)<10e-7) then
! uA=0.0D0
! end if
! if (abs(vA)<10e-7) then
! vA=0.0D0
! end if
! 
! end if
! 
! qA(1,j)=hA
! qA(2,j)=uA
! qA(3,j)=vA
! 
! END DO
! 
! 
! END SUBROUTINE OUTFLOW0_xi

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

SUBROUTINE OUTFLOWN_xi(fopt,tipo,MC,Nx,Ny,Fr2,dep,hs,Ns,t,dt,q,zt,ep_x,qA,zA)
!Subrutina que fija una altura de agua en el borde xi=N del dominio

USE coords
USE time0
!Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
implicit none
real (kind=8):: C,t,dt,Fr2,epA, hL, uL,vL, epL, epxL,epyL, RL,dep, &
zepL, epxA,epyA,zepA,tauL,Rmenos, hb,hA,zmin,uA,vA, Rmas, U1, alfa, &
uU,vU,U1u,hU,RmasU,uR,vR, Upar, UparR, UparL
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,v,zep_x,epx,epy
real (kind=8), dimension(Ns,2)::hs
real (kind=8), dimension(Nx,Ny)::MC
borde=2
!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(hs(:,1))) THEN
!      hb=hs(Ns,2)
     print*, 'Faltan datos hidrograma'
     stop
     !hb=hs(Ns,2)
ELSE
    call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hb)
    !print*,hb
!     print*,'Ns=',Ns
    !hb es la profundidad 
    !hs=(t,h)
END IF

epA=Nx
zmin=minval(zt(3,3:Ny+2))

DO j=1,Ny

zA(j)=zt(Nx+2,j+2)

if ((hb).le.zA(j)) then !nodo seco
hA=0.0D0
uA=0.0D0
vA=0.0D0

else

epx=ep_x(1,:,j)
epy=ep_x(2,:,j)

call interp1(Nx,coordxi,epx,epA,epxA)
call interp1(Nx,coordxi,epy,epA,epyA)

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO

! Caracteristica Saliente
!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
v=q(3,:,j)
if (h(Nx) == 0.0D0) then  !Caso seco
hL=0.0D0
uL=0.0D0
vL=0.0D0
RL=0.0D0
!Rmas=0.0D0
else
  !print*, j
  !call fzeroEP(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epL,epxL,epyL,uL,vL,hL,zepL)
  call fzeroEP_brent(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epL,epxL,epyL,uL,vL,hL,zepL)
  RL=(uL*epxL+vL*epyL)+2.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))
end if
     
if (fopt==0) then
tauL=0.0D0
Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
else
    if(hL/=0.0D0) then
    C=MC(Nx,j)
    call  tauU(tipo,C,Fr2,hL,uL,vL,tauL)
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+&
    zepA*(epxA**2.0D0+epyA**2.0D0))-tauL/hL*dt*(epxL+epyL)
    else
    tauL=0.0D0
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if
end if

!Situacion Borde
!Informacion: altura fija.
!Info en t=0

hU=q0_global(1,Nx,j)
uU=q0_global(2,Nx,j)
vU=q0_global(3,Nx,j)

U1u=uU*epxA+vU*epyA

RmasU=U1u+2.0D0*sqrt(hU/Fr2)*sqrt(epxA**2.0D0+epyA**2.0D0)

Rmenos=RmasU-4.0D0*sqrt(epxA**2.0D0+epyA**2.0D0)*sqrt(hb/Fr2)

hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)
U1=(Rmas+Rmenos)/2.0D0
alfa=angulo2(j)

uR=U1*cos(alfa)
vR=U1*sin(alfa)

UparL=-uL*sin(alfa)+vL*cos(alfa)
UparR=-uR*sin(alfa)+vR*cos(alfa)

if (U1>0.0D0) then
Upar=UparL
else
Upar=UparR
end if

uA=U1*cos(alfa)-Upar*sin(alfa)
vA=U1*sin(alfa)+Upar*cos(alfa)

if (abs(uA)<10e-7) then
uA=0.0D0
end if
if (abs(vA)<10e-7) then
vA=0.0D0
end if



end if
qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA

!print*, j,hA,uA,hb,Rmas,Rmenos, U1u

END DO

! pause


END SUBROUTINE OUTFLOWN_xi

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

SUBROUTINE OUTFLOWN_V(fopt,tipo,MC,Nx,Ny,Fr2,dep,vs,Ns,t,dt,q,zt,ep_x,qA,zA)
!Subrutina que fija una velocidad en el borde xi=N del dominio

USE coords
USE time0
!Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
implicit none
real (kind=8):: C,t,dt,Fr2,epA, hL, uL,vL, epL, epxL,epyL, RL,dep, &
zepL, epxA,epyA,zepA,tauL,Rmenos, hb,hA,zmin,uA,vA, Rmas, U1, alfa,&
ub, hU,uU,vU,U1u,RmasU, uR, vR, UparL, UparR, Upar
integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,v,zep_x,epx,epy
real (kind=8), dimension(Ns,2)::vs
real (kind=8), dimension(Nx,Ny)::MC
borde=2
!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(vs(:,1))) THEN
!      hb=hs(Ns,2)
     print*, 'Faltan datos hidrograma'
     stop
ELSE
    call interp1(Ns,vs(:,1),vs(:,2),(t+dt),ub)
    !hb es la cota (hs+zmin)
    !vs=(t,v)
END IF

epA=Nx+0.5D0
zmin=minval(zt(3,3:Ny+2))


DO j=1,Ny



zA(j)=zt(3,j+2)
hb=100000D0
if ((hb).le.zA(j)) then !nodo seco
hA=0.0D0
uA=0.0D0
vA=0.0D0

else

epx=ep_x(1,:,j)
epy=ep_x(2,:,j)

call interp1(Nx,coordxi,epx,epA,epxA)
call interp1(Nx,coordxi,epy,epA,epyA)

DO i=1,Nx
zep_x(i)=(zt(i+3,j+2)-zt(i+2,j+2))/dep	!j+2 porque incluye celdas ficticias
END DO

! Caracteristica Saliente
!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
v=q(3,:,j)
if (h(Nx) == 0.0D0) then  !Caso seco
hL=0.0D0
uL=0.0D0
vL=0.0D0
RL=0.0D0
!Rmas=0.0D0
else
  
  !print*, j
  
  call fzeroEP(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epL,epxL,epyL,uL,vL,hL,zepL)
  RL=(uL*epxL+vL*epyL)+2.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))
end if
     
if (fopt==0) then
tauL=0.0D0
Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
else
    if(hL/=0.0D0) then
    C=MC(Nx,j)
    call  tauU(tipo,C,Fr2,hL,uL,vL,tauL)
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+&
    zepA*(epxA**2.0D0+epyA**2.0D0))-tauL/hL*dt*(epxL+epyL)
    else
    tauL=0.0D0
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if
end if

!Situacion Borde
!Informacion saliente: velocidad fija U
vA=0.0D0
uA=ub
! Info en t=0
hU=q0_global(1,Nx,j)
uU=q0_global(2,Nx,j)
vU=q0_global(3,Nx,j)

U1u=uU*epxA+vU*epyA
RmasU=U1u+2.0D0*sqrt(hU/Fr2)*sqrt(epxA**2.0D0+epyA**2.0D0)

U1=uA*epxA+vA*epyA
Rmenos=2.0D0*U1-RmasU
!print*,U1

hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)
!U1=(Rmas+Rmenos)/2.0D0
alfa=angulo2(j)

uR=U1*cos(alfa)
vR=U1*sin(alfa)

UparL=-uL*sin(alfa)+vL*cos(alfa)
UparR=-uR*sin(alfa)+vR*cos(alfa)

if (U1>0.0D0) then
Upar=UparL
else
Upar=UparR
end if

uA=U1*cos(alfa)-Upar*sin(alfa)
vA=U1*sin(alfa)+Upar*cos(alfa)

if (abs(uA)<10e-7) then
uA=0.0D0
end if
if (abs(vA)<10e-7) then
vA=0.0D0
end if

! print*,U1,uA,ub, hA,hL
! pause

end if

qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA

!print*, hA, uA

END DO


END SUBROUTINE OUTFLOWN_V