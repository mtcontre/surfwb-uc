SUBROUTINE INFLOW0_eta(pasoRK,fopt,tipo,MC,Nx,Ny,Fr2,dep,Qsx,Qsy,etas3,timeS,Ns,t,dt,q,zt,ep_x,ep2_x,qA,zA)
!Basado en la rutina de Brett Sanders

!Subrutina que ingresa un q(t) unitario en cada nodo, en el borde xi=0, 
!Fijarse en el signo del caudal!!!
!Si q(t) entra qx o qy deben tener signo mas
!Si q(t) sale qx o qy deben apuntar hacia afuera, signo neg

USE coords
USE time0
!Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
implicit none
real (kind=8):: C,t,dt,Fr2,epA, hR, uR,vR, epR, epxR,epyR, RR,dep, zepR, &
  epxA,epyA,epxA2,epyA2,zepA,tauR,Rmenos,hA,zmin,uA,vA, Rmas,RmenosU,U1, &	
  alfa,a1,a2,a3,a4,K,qp,U1u,uU,vU,hU, JacA, cero, us, vs, hs
real (kind=8)::etaL,uL,vL,hL,RL, epL,epxL,epyL, zepL,C0, h0, tauL, qpar,uRper,vRpar,beta, pi,UparR,Upar, UparL, Uper 

integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde, pasoRK
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x, ep2_x
real (kind=8), dimension(3,Nx)::qA
real (kind=8), dimension(Nx)::zA,qx,qy, eta1
real (kind=8), dimension(Ny)::h,u,v,zep_x,epx,epy,epx2, epy2
real (kind=8), dimension(Nx,Ns)::Qsx,Qsy, etas3
real (kind=8), dimension(Ns)::timeS,qs1,qs2,etas
real (kind=8), dimension(Nx,Ny)::MC

real (kind=8):: zumbral


zumbral=-25.0D0+103.416330D0!Cambiar este valor con el Zmin de la batimetría
borde=3
!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(timeS)) THEN
!      hb=hs(Ns,2)
     print*, 'Faltan datos hidrograma'
     stop

ELSE

    !En el caso de que se tenga un qx y qy para cada nodo del borde
    !esto se debe hacer para cada j, usando
    Do j=1,Nx

    qs1=Qsx(j,:)  
    qs2=Qsy(j,:)
    etas=etas3(j,:)  
    alfa=angulo3(j)
    call interp1(Ns,timeS(:),qs1(:),(t+dt),qx(j))
    call interp1(Ns,timeS(:),qs2(:),(t+dt),qy(j))
    call interp1(Ns,timeS(:),etas(:),(t+dt),eta1(j))
    end Do

END IF

!Distribucion Q en q*
!Esto es por si no se tiene una distribucion en cada nodo

!Borde

epA=0.0D0
zmin=minval(zt(3:Nx+2,3))

DO j=1,Nx

alfa=angulo3(j)

zA(j)=zt(j+2,3)

h0=abs(zA(j)-103.416330D0)+1.4D0	!Cambiar considerando el Zmin y la marea máxima

etaL=eta1(j)

epx=ep_x(1,j,:) !Eta
epy=ep_x(2,j,:) !Eta

epx2=ep2_x(1,j,:) !Xi
epy2=ep2_x(2,j,:) !Xi

! print*,'j=',j
!print*,'epx2=',epx2

call interp1(Ny,coordeta,epx,epA,epxA) !Eta_x
call interp1(Ny,coordeta,epy,epA,epyA) !Eta_y

call interp1(Ny,coordeta,epx2,epA,epxA2) !Xi_x
call interp1(Ny,coordeta,epy2,epA,epyA2) !Xi_y

JacA=epxA2*epyA-epyA2*epxA
! ! print*,'xi_x=',epxA2,'xi_y=',epyA2,'eta_x=',epxA,'eta_y=',epyA
! print*,'JacA=',JacA
! ! print*,j
! pause


DO i=1,Ny
zep_x(i)=(zt(j+2,i+3)-zt(j+2,i+2))/dep	!j+2 porque incluye celdas ficticias
END DO

zepA=(zt(j+2,3)-zt(j+2,2))/dep 
!zepA=(zt(j+2,3)-zt(j+2,2))/dep 

! Caracteristica Saliente
!Rmenos
h=q(1,j,:)	
u=q(2,j,:)
v=q(3,j,:)


!Esto es para la modelacion de tsunami solamente!!
if (zA(j)>=zumbral) then !Borde abierto

qA(1,j)=h(1)
qA(2,j)=u(1)
qA(3,j)=v(1)
! print*, 'BA', j,zA(j)-104.01078D0
! pause

else !genabs


if (h(1) == 0.0D0) then  !Caso seco
hR=0.0D0
uR=0.0D0
vR=0.0D0
RR=0.0D0

else
  !print*, h(1),j
  call fzeroEP(j,borde,Fr2,Ny,dep,epx,epy,h,u,v,zep_x,dt,epR,epxR,epyR,uR,vR,hR,zepR)
  !call fzeroEP_brent(borde,Fr2,Ny,dep,epx,epy,h,u,v,zep_x,dt,epR,epxR,epyR,uR,vR,hR,zepR)

  ! Usando la direccion Eta como la perpendicular al borde
! print*, ZepR
! pause
  RR=(uR*epxR+vR*epyR)-2.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))
  
!   ! En la direccion normal al borde
!    uRper=uR*cos(alfa)+vR*sin(alfa)
!    vRpar=-uR*sin(alfa)+vR*cos(alfa)
!    
!    RR=uRper-2.0D0*sqrt(hR/Fr2)



!print*,'RR=', epR,RR, epxR, epyR, uR, vR, hR, 'UR=',uR*epxR+vR*epyR, 'Fr2=',Fr2
   !pause
end if
     
if (fopt==0) then
tauR=0.0D0
Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
else
    if(hR/=0.0D0) then
    C=MC(j,1)
    call  tauU(tipo,C,Fr2,hR,uR,vR,tauR)
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))-tauR/hR*dt*(epxR+epyR)
    else
    tauR=0.0D0
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if

end if

!Caracteristica Entrante



! En la perpendicular al borde
!qp=qy(j)*cos(alfa)+qx(j)*sin(alfa)

! En la direccion eta
qp=qx(j)*epxA+qy(j)*epyA
qpar=qx(j)*epxA2+qy(j)*epyA2

!---------------------------------------------
! Calculando RL con ambos datos desde afuera
!qpar=-qx(j)*sin(alfa)+qy(j)*cos(alfa)

 C0=sqrt(h0/Fr2)
epL=epA-dt*C0
hL=etaL+h0
if (hL<=0.0D0) then
hL=0.0D0
uL=0.0D0
vL=0.0D0
else
uL=qx(j)/hL
vL=qy(j)/hL
end if

call interp1(Ny,coordeta,zep_x,epL,zepL)
call interp1(Ny,coordeta,epx,epL,epxL)
call interp1(Ny,coordeta,epy,epL,epyL)


! Usando como perpendicular la direccion eta
RL=(uL*epxL+vL*epyL)+2.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))

! print*, RL
! pause

!Rmas
if (fopt==0) then
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(j,1)
  call tauU(tipo,C,Fr2,hL,uL,vL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if
Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*(epxL+epyL)*dt



!print*,'hL=',hL,'uL=',uL,'vL=',vL,'etaL=',eta1(j),'zA=',zA(j)-105.37483D0
! print*,'Rmenos=',Rmenos,'Rmas=',Rmas
!!--------------------------------------
!Metodologia Brett Sanders Alturas Fijas
!!--------------------------------------
! ! Completo con hR, uR, vR el q0 en la primera iteracion
! call q0(borde,pasoRK,hR,uR,vR,j)
! 
! !Solucion Ec. Orden 3
! hU=q0_global(1,j,1)
! uU=q0_global(2,j,1)
! vU=q0_global(3,j,1)
!  
! !U perpendicular U1 asumiendo que Uper=U1 en la direccion Xi
! U1u=uU*epxA+vU*epyA 
! RmenosU=U1u-2.0D0*sqrt(hU/Fr2)*sqrt(epxR**2.0D0+epyR**2.0D0)
! 
! RL=RmenosU+4.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))
! 
! if (fopt==0) then
! tauL=0.0D0
! else
!   if (hL/=0.0D0) then
!   C=MC(1,j)
!   call tau(tipo,C,Fr2,hL,uL,vL,tauL)
!   tauL=tauL/hL
!   else
!   tauL=0.0D0
!   end if
! end if
! Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*(epxL+epyL)*dt

! ! ---------------------------------------
! ! Usando metodología de Brett Sanders para caudales
! ! --------------------
! ! Completo con hR, uR, vR el q0 en la primera iteracion
! call q0(borde,pasoRK,hR,uR,vR,j)
! !----------
! !Solucion Ec. Orden 3
! hU=q0_global(1,j,1)
! uU=q0_global(2,j,1)
! vU=q0_global(3,j,1)
! 
! !U perpendicular U1 asumiendo que Uper=U1 en la direccion eta
! U1u=uU*epxA+vU*epyA
! RmenosU=U1u-2.0D0*sqrt((hU/Fr2)*(epxA**2.0D0+epyA**2.0D0))
! !print*,'U1u=',U1u, 'hU=',hU, 'epyA=',epyA,'epyR=',epyR
! 
! !En la perpendicular al borde
! ! U1u=uU*cos(alfa)+vU*sin(alfa)
! ! RmenosU=U1u-2.0D0*sqrt(h	U/Fr2)
! K=Fr2/(32.0D0*(epxA**2.0D0+epyA**2.0D0))
! a1=K
! a2=-K*RmenosU
! a3=-K*RmenosU**2.0D0
! if (qp>=0.0D0) then
! a4=K*RmenosU**3.0D0-qp
! else
! a4=K*RmenosU**3.0D0+qpabout:startpage
! end if
! 
! !a4=K*RmenosU**3.0D0-qp
! call grado3(a1,a2,a3,a4,RL)
! !Verificacion
!  cero=a1*RL**3.0D0+a2*RL**2.0D0+a3*RL+a4
!  
! print*,'0=', cero
! !Rmas
! hL=etaL+h0
! uL=qx(j)/hL
! vL=qy(j)/hL
! 
! if (fopt==0) then
! tauL=0.0D0about:startpage
! else
!   if (hL/=0.0D0) then
!   C=MC(1,j)
!   call tau(tipo,C,Fr2,hL,uL,vL,tauL)
!   tauL=tauL/hL
!   else
!   tauL=0.0D0
!   end if
! end if
! Rmas=RL-dt/Fr2*0.5D0*(zepL	*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*(epxL+epyL)*dt
! 


!!!---------------------------------------------------------
! Solucion:

!En la direccion perpendicular dominio curvilíneo Xi
hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)

! En la direccion perpendicular dominio real
!hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0
Uper=(Rmas+Rmenos)/2.0D0

! print*,'Rmas=', Rmas, 'Rmenos=',Rmenos!, 'RmenosU=',RmenosU
! print*, 'Rmas-Rmenos=',Rmas-Rmenos


UparL=uL*epxA2+vL*epyA2  !U1=u*xi_x+v*xi_y
UparR=uR*epxA2+vR*epyA2




if (Uper>0.0D0) then
Upar=UparL
else
Upar=UparR
end if

! if (qp<=0.0D0) then
uA=Upar*epyA/JacA+Uper*-epyA2/JacA !U1*eta_y/J-U2*xi_y/J * Uper=U2, Upar=U1
vA=Upar*-epxA/JacA+Uper*epxA2/JacA !U1*-eta_x/J+U2*xi_x/J



! else
! uA=Uper*-epyA2/JacA+Upar*epyA/JacA
! vA=Uper*epxA2/JacA+Upar*-epxA/JacA
! end if

! En la perpendicular al dominio
! vA=U1*cos(alfa)-Upar*sin(alfa)
! uA=U1*sin(alfa)-Upar*cos(alfa)

if (hA<1e-7) then
hA=0.0D0
uA=0.0D0
vA=0.0D0
qA(1,j)=0.0D0
qA(2,j)=0.0D0
qA(3,j)=0.0D0
else
qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA
end if


end if

! print*,'j=',j
! print*, 'hA=',hA, 'uA=',uA, 'vA=',vA
! print*, 'qxA=',hA*uA, 'qyA=',hA*vA, 'dh0=',hA-h0,'dhR=',hA-hR
! print*, 'qx=',qx(j), 'qy=',qy(j)!, 'uL=',uL,'vL=',vL
! ! print*,'hL=',hL,'hR=',hR
! ! print*, 'Upar=',Upar,'Uper=',Uper, 'JacA=',JacA
! print*, 'xi_x=',epxA2,'xi_y=',epyA2,'eta_x=',epxA,'eta_y=',epyA
! pause

END DO



END SUBROUTINE INFLOW0_eta


SUBROUTINE INFLOWN_eta(pasoRK,fopt,tipo,MC,Nx,Ny,Fr2,dep,Qsx,Qsy,etas4,timeS,Ns,t,dt,q,zt,ep_x,ep2_x,qA,zA)
!Basado en la rutina de Brett Sanders

!Subrutina que ingresa un q(t) unitario en cada nodo, en el borde xi=N, 
!Fijarse en el signo del caudal!!!
!Si q(t) entra qx o qy deben tener signo neg
!Si q(t) sale qx o qy deben apuntar hacia afuera, signo positivo

!Por el momento está hecha para recibir un Qtotal en la seccion
USE coords
USE time0
!Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
implicit none
real (kind=8):: C,t,dt,Fr2,epA, hL, uL,vL, epL, epxL,epyL, RL,dep, zepL, &
  epxA,epyA,epxA2,epyA2,zepA,tauL,Rmenos, qb,hA,zmin,uA,vA, Rmas,RmasU,U1,&
    alfa, AT,QT,uT,a1,a2,a3,a4,K,qp,U1u,uU,vU,hU, Upar,UparL, UparR, JacA, cero

real (kind=8)::etaR,uR,vR,hR,RR, epR,epxR,epyR, zepR,C0, h0, tauR, beta, pi,qpar,uLper,vLpar, Uper

integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde, pasoRK
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x, ep2_x
real (kind=8), dimension(3,Nx)::qA
real (kind=8), dimension(Nx)::zA, qx, qy, eta1
real (kind=8), dimension(Ny)::h,u,v,zep_x,epx,epy, epx2, epy2
real (kind=8), dimension(Ns,2)::qs
real (kind=8), dimension(Nx,Ns)::Qsx,Qsy, etas4
real (kind=8), dimension(Ns)::timeS,qs1,qs2,etas
real (kind=8), dimension(Nx,Ny)::MC

real (kind=8):: zumbral
zumbral=-20.0D0+103.416330D0!Cambiar
borde=4
!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(timeS)) THEN
!      hb=hs(Ns,2)
     print*, 'Faltan datos hidrograma'
     stop
ELSE
!     call interp1(Ns,qs(:,1),qs(:,2),(t+dt),qb)
!     qb es el Qtotal entrando en la direccion -X
    
    !En el caso de que se tenga un qx y qy para cada nodo del borde
    !esto se debe hacer para cada j, usando
    Do j=1,Nx
    qs1=Qsx(j,:) !El de direccion y es el perpendicular al borde
    qs2=Qsy(j,:)
    etas=etas4(j,:)
    alfa=angulo4(j)
    call interp1(Ns,timeS(:),qs1(:),(t+dt),qx(j))
    call interp1(Ns,timeS(:),qs2(:),(t+dt),qy(j))
    call interp1(Ns,timeS(:),etas(:),(t+dt),eta1(j))
    end Do

END IF



!Borde

epA=Ny
zmin=minval(zt(3:Nx+2,Ny+2))

DO j=1,Nx

alfa=angulo4(j)

zA(j)=zt(j+2,Ny+2)
h0=abs(zA(j)-103.416330D0)+1.4D0!Cambiar en cada caso de batimetría y marea
etaR=eta1(j)


! print*,h0,etaR
! pause

epx=ep_x(1,j,:)
epy=ep_x(2,j,:)
epx2=ep2_x(1,j,:) !Xi
epy2=ep2_x(2,j,:) !Xi

call interp1(Ny,coordeta,epx,epA,epxA)
call interp1(Ny,coordeta,epy,epA,epyA)
call interp1(Ny,coordeta,epx2,epA,epxA2) !Xi_x
call interp1(Ny,coordeta,epy2,epA,epyA2) !Xi_y

JacA=epxA2*epyA-epyA2*epxA

DO i=1,Ny
zep_x(i)=(zt(j+2,i+3)-zt(j+2,i+2))/dep	!j+2 porque incluye celdas ficticias
END DO

zepA=(zt(j+2,Ny+3)-zt(j+2,Ny+2))/dep
! print*, 'zepA=',zepA
! Caracteristica Saliente



!Rmenos
h=q(1,j,:)	
u=q(2,j,:)
v=q(3,j,:)

if (zA(j)>=zumbral) then !tsunami !borde abierto
!Esto es para la modelacion de tsunami solamente!!
!if (zA(j)>=zumbral) then
qA(1,j)=h(Ny-1)
qA(2,j)=u(Ny-1)
qA(3,j)=v(Ny-1)

else !genabs


if (h(Ny) == 0.0D0) then  !Caso seco
hL=0.0D0
uL=0.0D0
vL=0.0D0
RL=0.0D0
!Rmas=0.0D0
else
  !print*, j
  call fzeroEP(j,borde,Fr2,Ny,dep,epx,epy,h,u,v,zep_x,dt,epL,epxL,epyL,uL,vL,hL,zepL)
  !call fzeroEP_brent(borde,Fr2,Ny,dep,epx,epy,h,u,v,zep_x,dt,epL,epxL,epyL,uL,vL,hL,zepL)

  ! Usando la direccion Eta como la perpendicular al borde
  RL=(uL*epxL+vL*epyL)+2.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))

!   ! En la direccion normal al borde
!    uLper=uL*cos(alfa)+vL*sin(alfa)
!    vLpar=-uL*sin(alfa)+vL*cos(alfa)
!    
!    RL=uLper-2.0D0*sqrt(hL/Fr2)



!   print*,'RL=', epL,RL, epxL, epyL, uL, vL, hL, 'UL=',uL*epxL+vL*epyL, 'Fr2=',Fr2
!   pause
end if
     
if (fopt==0) then
tauL=0.0D0
Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
else
    if(hL/=0.0D0) then
    C=MC(j,Ny)
    call  tauU(tipo,C,Fr2,hL,uL,vL,tauL)
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))-tauL/hL*dt*(epxL+epyL)
    else
    tauL=0.0D0
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if
end if


!Caracteristica Entrante

! En la perpendicular al borde
!qp=qy(j)*cos(alfa)+qx(j)*sin(alfa)

! En la direccion eta
qp=qx(j)*epxA+qy(j)*epyA
qpar=qx(j)*epxA+qy(j)*epyA
!---------------------------------------------
! ! Calculando RR con ambos datos desde afuera

 C0=sqrt(h0/Fr2)
epR=epA+dt*C0
hR=etaR+h0

if (hR<=0.0D0) then
hR=0.0D0
uR=0.0D0
vR=0.0D0
else
uR=qx(j)/hR
vR=qy(j)/hR
end if

call interp1(Ny,coordeta,zep_x,epR,zepR)
call interp1(Ny,coordeta,epx,epR,epxR)
call interp1(Ny,coordeta,epy,epR,epyR)

! print*, 'epxR=',epxR,'epyR=',epyR
RR=uR*epxR+vR*epyR-2.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))

!Rmenos
if (fopt==0) then
tauR=0.0D0
else
  if (hR/=0.0D0) then
  C=MC(j,Ny)
  call tauU(tipo,C,Fr2,hR,uR,vR,tauR)
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if
Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*(epxR+epyR)*dt

! print*,'hR=',hR,'uR=',uR,'vR=',vR,'etaR=',eta1(j),'zA=',zA(j)

! ! ---------------------------------------
! ! Usando metodología de Brett Sanders
! ! --------------------
!!--------------------------------------
!Metodologia Brett Sanders Alturas Fijas
!!--------------------------------------
! ! Completo con hR, uR, vR el q0 en la primera iteracion
! call q0(borde,pasoRK,hL,uL,vL,j)
! 
! !Solucion Ec. Orden 3
! hU=q0_global(1,j,Ny)
! uU=q0_global(2,j,Ny)
! vU=q0_global(3,j,Ny)
!  
! !U perpendicular U1 asumiendo que Uper=U1 en la direccion Xi
! U1u=uU*epxA+vU*epyA 
! RmasU=U1u+2.0D0*sqrt(hU/Fr2)*sqrt(epxL**2.0D0+epyL**2.0D0)
! 
! RR=RmasU-4.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))
! 
! !Rmenos
! if (fopt==0) then
! tauR=0.0D0
! else
!   if (hR/=0.0D0) then
!   C=MC(1,j)
!   call tau(tipo,C,Fr2,hR,uR,vR,tauR)
!   tauR=tauR/hR
!   else
!   tauR=0.0D0
!   end if
! end if
! Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*(epxR+epyR)*dt

! print*,'hR=',hR,'uR=',uR,'vR=',vR,'etaR=',eta1(j),'zA=',zA(j)

!----------------------
! Metodologia Caudales

! ! Completo con hR, uR, vR el q0 en la primera iteracion
! call q0(borde,pasoRK,hL,uL,vL,j)
! !-----------------
! !Solucion Ec. Orden 3
! hU=q0_global(1,j,Ny)
! uU=q0_global(2,j,Ny)
! vU=q0_global(3,j,Ny)
! 
! !U perpendicular U1 asumiendo que Uper=U1 en la direccion eta
! U1u=uU*epxA+vU*epyA
! RmasU=U1u+2.0D0*sqrt((hU/Fr2)*(epxA**2.0D0+epyA**2.0D0))
! ! print*,'0=', cero, qp
! ! print*,'U1u=',U1u, 'hU=',hU, 'epyA=',epyA,'epyL=',epyL
! 
! !En la perpendicular al borde
! ! U1u=uU*cos(alfa)+vU*sin(alfa)
! ! RmenosU=U1u-2.0D0*sqrt(hU/Fr2)
! 
! K=Fr2/(32.0D0*(epxA**2.0D0+epyA**2.0D0))
! a1=K
! a2=-K*RmasU
! a3=-K*RmasU**2.0D0
! if (qp<=0.0D0) then
! a4=K*RmasU**3.0D0-qp
! else
! a4=K*RmasU**3.0D0+qp
! end if
! 
! call grado3(a1,a2,a3,a4,RR)
! !Verificacion
!  cero=a1*RR**3.0D0+a2*RR**2.0D0+a3*RR+a4
! ! print*,'0=', cero, qp
! 
! !Rmenos
! hR=etaR+h0
! uR=qx(j)/hR
! vR=qy(j)/hR
! 
! if (fopt==0) then
! tauR=0.0D0
! else
!   if (hR/=0.0D0) then
!   C=MC(1,j)
!   call tau(tipo,C,Fr2,hR,uR,vR,tauR)
!   tauR=tauR/hR
!   else
!   tauR=0.0D0
!   end if
! end if
! Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*(epxR+epyR)*dt
!!!!
!!----------------------------------------------------------
!Solucion

!En la direccion perpendicular dominio curvilíneo Xi
hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)

! En la direccion perpendicular dominio real
!hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0

Uper=(Rmas+Rmenos)/2.0D0
! print*,'Rmas=', Rmas, 'Rmenos=',Rmenos, 'Uper=',Uper

! uR=qx(j)/hA
! vR=qy(j)/hA
!Uper=U1

UparL=uL*epxA2+vL*epyA2  !U1=u*xi_x+v*xi_y
UparR=uR*epxA2+vR*epyA2

!UparL=-uL*cos(alfa)-vL*sin(alfa)
!UparR=-uR*cos(alfa)-vR*sin(alfa)

! print*, 'Uper=', Uper !, 'UparL=',UparL,'UparR=',UparR,'Upar=',Upar 
if (Uper>0.0D0) then
Upar=UparR
! print*, 'Upar=', Upat
else
Upar=UparL 
! 
end if
! print*, 'Upar=', Upar

!if (qp>=0.0D0) then 
uA=Upar*epyA/JacA+Uper*-epyA2/JacA !-U1*epyA2/JacA+Upar*epyA/JacA
vA=Upar*-epxA/JacA+Uper*epxA2/JacA


!else
! uA=Upar*epyA/JacA+Uper*-epyA2/JacA
! vA=Upar*-epxA/JacA+Uper*epxA2/JacA
! end if

! print*,epxA2, epxA


! print*,'U1=',U1, 'Upar=',Upar!,'R=', UparR, 'L=',UparL !, 'RmenosU=',RmenosU
! pause

! En la perpendicular al dominio
! vA=-U1*cos(alfa)+Upar*sin(alfa)
! uA=U1*sin(alfa)-Upar*cos(alfa)
! vA=-U1*cos(beta)-Upar*sin(beta)
! uA=-U1*sin(beta)+Upar*cos(beta)


if (hA<1e-7) then
uA=0.0D0
vA=0.0D0
hA=0.0D0
qA(1,j)=0.0D0
qA(2,j)=0.0D0
qA(3,j)=0.0D0
else
qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA
end if


end if

!print*, hU,uU,j

! print*,j 
! print*, 'hA=',hA, 'uA=',uA, 'vA=',vA
! print*, 'qxA=',hA*uA, 'qyA=',hA*vA, 'dh0=',h0-hA,'dhL=',hA-hL
! print*, 'qx=',qx(j),'qy=',qy(j), 'etaR=',etaR
! pause

END DO




END SUBROUTINE INFLOWN_eta
