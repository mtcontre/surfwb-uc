SUBROUTINE INFLOW0_xi(pasoRK,fopt,tipo,MC,Nx,Ny,Fr2,dep,Qsx,Qsy,etas1,timeS,Ns,t,dt,q,zt,ep_x,ep2_x,qA,zA)
!Basado en la rutina de Brett Sanders

!Subrutina que ingresa un q(t) unitario en cada nodo, en el borde xi=0, 
!Fijarse en el signo del caudal!!!
!Si q(t) entra qx o qy deben tener signo mas
!Si q(t) sale qx o qy deben apuntar hacia afuera, signo neg

USE coords
USE time0
!Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0
implicit none
real (kind=8):: C,t,dt,Fr2,epA, hR, uR,vR, epR, epxR,epyR, RR,dep, zepR, epxA,epyA,epxA2,epyA2,zepA,tauR,Rmenos,hA,zmin,uA,vA, Rmas,RmenosU,U1, alfa,a1,a2,a3,a4,K,qp,U1u,uU,vU,hU, cero, JacA, us, vs, hs
real (kind=8)::etaL,uL,vL,hL,RL, epL,epxL,epyL, zepL,C0, h0, tauL, qpar, Uper,Upar, UparL, UparR 

integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde, pasoRK
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x, ep2_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA,qx,qy, eta1
real (kind=8), dimension(Nx)::h,u,v,zep_x,epx,epy, epx2, epy2
real (kind=8), dimension(Ny,Ns)::Qsx,Qsy,etas1
real (kind=8), dimension(Ns)::timeS,qs1,qs2,etas
real (kind=8), dimension(Nx,Ny)::MC
real (kind=8):: zumbral, x1,x2, x3, qx2, qy2, u0,v0

borde=1
zumbral=0.0D0+103.416330D0 !Cambiar este valor
!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(timeS)) THEN
!      hb=hs(Ns,2)
!      print*, 'Faltan datos hidrograma'
!      stop
    Do j=1,Ny
    qx(j)=0.0D0
    qy(j)=0.0D0
    eta1(j)=0.0D0
    end do
    

ELSE

    !En el caso de que se tenga un qx y qy para cada nodo del borde
    !esto se debe hacer para cada j, usando
    Do j=1,Ny
    qs1=Qsx(j,:)
    qs2=Qsy(j,:)
    etas=etas1(j,:)
    alfa=angulo1(j)
!     print*, qs1
!     pause
!     print*, etas
!     pause
!     print*, timeS
! pause
    call interp1(Ns,timeS(:),qs1(:),(t+dt),qx(j))
    call interp1(Ns,timeS(:),qs2(:),(t+dt),qy(j))
    call interp1(Ns,timeS(:),etas(:),(t+dt),eta1(j))
!     print*,t+dt*cos(alfa),eta1(j)
  ! print*, eta1(j)
!     pause
    ! Leer h, u y v 
    end Do
!     print*, 't=t',t
!      print*, 'etas(:)=',etas1(1,:)
   !   pause
END IF

!Borde

epA=0.0D0

zmin=minval(zt(3,3:Ny+2))

DO j=1,Ny

!Borde
zA(j)=zt(3,j+2)
h0=abs(zA(j)-103.416330D0)+1.4D0!Cambiar este valor sumar la marea alta
etaL=eta1(j)
! print*,h0,etaL
! pause
alfa=angulo1(j)

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
  zepR=0.0D0
!   ! Usando la direccion Xi como la perpendicular al borde
    RR=(uR*epxR+vR*epyR)-2.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))
  
!   ! En la direccion normal al borde
!    uRper=uR*cos(alfa)+vR*sin(alfa)
!    vRpar=-uR*sin(alfa)+vR*cos(alfa)
!    
!    RR=uRper-2.0D0*sqrt(hR/Fr2)
! print*,j
!      print*,'RR=', epR,RR, epxR, epyR, uR, vR, hR, 'UR=',uR*epxR+vR*epyR, 'zepR=',zepR,'zepA=',zepA
!      pause
     
  
end if
     
if (fopt==0) then
tauR=0.0D0
Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
else
    if(hR/=0.0D0) then
    C=MC(1,j)
    call  tauU(tipo,C,Fr2,hR,uR,vR,tauR)
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))-tauR/hR*dt*(epxR+epyR)
    else
    tauR=0.0D0
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if

end if

! Caracteristica Entrante

!----------------------------------
! En la perpendicular al borde
!qp=qx(j)*cos(alfa)+qy(j)*sin(alfa)


! En la direccion Xi
qp=qx(j)*epxA+qy(j)*epyA
qpar=qx(j)*epxA2+qy(j)*epyA2

! ! ---------------------------------------------
! Calculando RL con ambos datos desde afuera
!qpar=-qx(j)*sin(alfa)+qy(j)*cos(alfa)


 C0=sqrt(h0/Fr2)
epL=epA-dt*C0
hL=etaL+h0
! print*, 'L=',zA(j)-207.9032D0, 'h0=',h0, etaL, hL
! pause
if (hL<=0.0D0) then
hL=0.0D0
uL=0.0D0
vL=0.0D0
else
uL=qx(j)/hL  
vL=qy(j)/hL  

end if

call interp1(Nx,coordxi,zep_x,epL,zepL)
call interp1(Nx,coordxi,epx,epL,epxL)
call interp1(Nx,coordxi,epy,epL,epyL)
 zepL=0.0D0
! Usando como perpendicular la direccion Xi
RL=(uL*epxL+vL*epyL)+2.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))

!Usando la recta perpendicular al borde
!RL=uL+2.0D0*sqrt(hL/Fr2)

!Rmas
if (fopt==0) then
tauL=0.0D0
else
  if (hL/=0.0D0) then
  C=MC(1,j)
  call tauU(tipo,C,Fr2,hL,uL,vL,tauL)
  tauL=tauL/hL
  else
  tauL=0.0D0
  end if
end if
Rmas=RL-dt/Fr2*0.5D0*(zepL*epxL**2.0D0+zepA*epxA**2.0D0)-tauL*(epxL+epyL)*dt

! print*,'hL=',hL,'uL=',uL,'vL=',vL!,'etaL=',eta1(j),'zA=',zA(j)
! print*, 'Rmas=',Rmas,'RL=',RL!, 'epxL=',epxL,'epyL=',epyL
!  print*, 'zepL=',zepL,'zepA=',zepA, 'tauL=',tauL, 'epL=',epL
!  pause

! ! ---------------------------------------
! ! Usando metodología de Brett Sanders
! ! ---------------------------------------
! Metodologia Alturas
! -----------------------------------------
! ! Completo con hR, uR, vR el q0 en la primera iteracion
! call q0(borde,pasoRK,hR,uR,vR,j)
! 
! !Solucion Ec. Orden 3
! hU=q0_global(1,1,j)
! uU=q0_global(2,1,j)
! vU=q0_global(3,1,j)
!  
! !U perpendicular U1 asumiendo que Uper=U1 en la direccion Xi
! U1u=uU*epxA+vU*epyA 
! RmenosU=U1u-2.0D0*sqrt(hU/Fr2)*sqrt(epxA**2.0D0+epyA**2.0D0)
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

! Metodolodia Caudales
! --------------------
! ! Completo con hR, uR, vR el q0 en la primera iteracion
! call q0(borde,pasoRK,hR,uR,vR,j)
! 
! !Solucion Ec. Orden 3
! hU=q0_global(1,1,j)
! uU=q0_global(2,1,j)
! vU=q0_global(3,1,j)
! 
! !U perpendicular U1 asumiendo que Uper=U1 en la direccion Xi
! U1u=uU*epxA+vU*epyA 
! RmenosU=U1u-2.0D0*sqrt(hU/Fr2)*sqrt(epxA**2.0D0+epyA**2.0D0)
! 
! !En la perpendicular al borde
! ! U1u=uU*cos(alfa)+vU*sin(alfa)
! ! RmenosU=U1u-2.0D0*sqrt(hU/Fr2)
! 
! 
! !qp=qx(j)*cos(alfa)+qy(j)*sin(alfa)
! ! print*, 'qp=',qp,qx(j),qy(j)
! K=Fr2/(32.0D0*(epxA**2.0D0+epyA**2.0D0))
! a1=K
! a2=-K*RmenosU
! a3=-K*RmenosU**2.0D0
! if (qp>=0.0D0) then
! a4=K*RmenosU**3.0D0-qp
! else
! a4=K*RmenosU**3.0D0+qp
! end if
! !a4=K*RmenosU**3.0D0-qp
! call grado3(a1,a2,a3,a4,RL)
! !Verificacion
!  cero=a1*RL**3.0D0+a2*RL**2.0D0+a3*RL+a4
! ! print*,'0=', cero, qp
! 
! !Rmas
! hL=etaL+h0
! uL=qx(j)/hL
! vL=qy(j)/hL
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
! ! ! 
!!!---------------------------------------------------------
! Solucion:
!En la direccion perpendicular dominio curvilíneo Xi
hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)
! En la direccion perpendicular dominio real
!hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0
Uper=(Rmas+Rmenos)/2.0D0
! print*,'Rmas=', Rmas, 'Rmenos=',Rmenos!, 'hA=',hA,'U1A=',U1
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

! En la direccion Xi
! if (qp<=0.0D0) then
! uA=Uper*epyA2/JacA+Upar*-epyA/JacA !U1*eta_y/J-Upar*Xiy/Jac
! vA=Uper*-epxA2/JacA+Upar*epxA/JacA
! else
uA=Uper*epyA2/JacA+Upar*-epyA/JacA !U1*eta_y/J-U2*Xiy/J
vA=Uper*-epxA2/JacA+Upar*epxA/JacA !U1*-eta_x/J+U2*xi_x/J
! end if


! print*, 'U1=', U1, 'UparL=',UparL,'UparR=',UparR,'Upar=',Upar 

! En la perpendicular al dominio
!uA=U1*cos(alfa)-Upar*sin(alfa)
!vA=U1*sin(alfa)+Upar*cos(alfa)


if (hA<1e-7) then
qA(1,j)=0.0D0
qA(2,j)=0.0D0
qA(3,j)=0.0D0
else
qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA
end if

end if

! ! ! ! ! !

! print*,'hL=',hL,'hR=',hR
! print*, 'hA=',hA, 'uA=',uA, 'vA=',vA
! print*, 'qxA=',hA*uA, 'qyA=',hA*vA!, 'dh0=',hA-h0,'dhR=',hA-hR
! print*, 'qx=',qx(j), 'qy=',qy(j)!, 'etaL=',etaL
! 
! print*, 'dh0=',hA-h0,'dhR=',hA-hR, 'dhL=',hL-hA
! ! print*, 'etaL=',etaL
! if (j==1) then
! print*,j
! ! print*,'Rmas=', Rmas, 'Rmenos=',Rmenos
! ! print*, 'Rmas-Rmenos=',Rmas-Rmenos
! ! print*, 'Upar=',Upar,'Uper=',Uper
! ! print*,'hL=',hL,'hR=',hR
! print*, 'hA=',hA, 'uA=',uA, 'vA=',vA
! print*, 'qxA=',hA*uA, 'qyA=',hA*vA, 'dh0=',hA-h0,'dhR=',hA-hR
! print*, 'xi_x=',epxA,'xi_y=',epyA,'eta_x=',epxA2,'eta_y=',epyA2
! print*, 'qx=',qx(j), 'qy=',qy(j)!, 'etaL=',etaL
! print*,'JacA=',JacA
! ! print*,j
! 
! ! print*, 'dh0=',hA-h0,'dhR=',hA-hR, 'dhL=',hL-hA
! ! print*, 'etaL=',etaL
! pause
! end if


END DO
!pause
! print*, 'hA=',hA, 'uA=',uA, 'vA=',vA
! print*, 'qxA=',hA*uA, 'qxR=',hR*uR, 'dh=',hA-hU 
! pause
! print*, 'hA=',hA, 'uA=',uA, 'vA=',vA
! print*, 'qA=',qA(1,Ny),qA(2,Ny),qA(3,Ny)
! print*, 'hR=',hR, 'uR=',uR !hA*uA, 'qxR=',hR*uR, 'dh=',hA-hU 
! pause



END SUBROUTINE INFLOW0_xi

!---------------------------------------------------------------------------

SUBROUTINE INFLOWN_xi(pasoRK,fopt,tipo,MC,Nx,Ny,Fr2,dep,Qsx,Qsy,etas2,timeS,Ns,t,dt,q,zt,ep_x,ep2_x,qA,zA)
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
real (kind=8):: C,t,dt,Fr2,epA, hL, uL,vL, epL, epxL,epyL, RL,dep, zepL, epxA,epyA,epxA2, epyA2,zepA,tauL,Rmenos, qb,hA,zmin,uA,vA, Rmas,RmasU,U1, alfa, AT,QT,uT,a1,a2,a3,a4,K,qp,U1u,uU,vU,hU,beta, JacA

real (kind=8):: etaR, hR, uR, vR, RR, epR, epxR, epyR, zepR, C0,h0, tauR, qpar, Uper, Upar, UparL, UparR 

integer:: Nx,Ny,j,i, Ns, fopt, tipo, borde, pasoRK
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x,ep2_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA, qx, qy, eta2
real (kind=8), dimension(Nx)::h,u,v,zep_x,epx,epy, epx2, epy2
real (kind=8), dimension(Ns,2)::qs
real (kind=8), dimension(Ny,Ns)::Qsx,Qsy, etas2
real (kind=8), dimension(Ns)::timeS,qs1,qs2, etas
real (kind=8), dimension(Nx,Ny)::MC

real (kind=8):: zumbral
zumbral=0.0D0+105.37483D0

borde=2
!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(timeS)) THEN
!      hb=hs(Ns,2)
     print*, 'dt=',dt,t,'Faltan datos hidrograma'
     stop
ELSE
!     call interp1(Ns,qs(:,1),qs(:,2),(t+dt),qb)
!     qb es el Qtotal entrando en la direccion -X
    
    !En el caso de que se tenga un qx y qy para cada nodo del borde
    !esto se debe hacer para cada j, usando
    Do j=1,Ny
    qs1=Qsx(j,:)
    qs2=Qsy(j,:)
    etas=etas2(j,:)
    call interp1(Ns,timeS(:),qs1(:),(t+dt),qx(j))
    call interp1(Ns,timeS(:),qs2(:),(t+dt),qy(j))
    call interp1(Ns,timeS(:),etas(:),(t+dt),eta2(j))

    end Do

END IF

!Borde

epA=Nx

zmin=minval(zt(Nx+2,3:Ny+2))

DO j=1,Ny

!Borde
zA(j)=zt(Nx+2,j+2)
h0=abs(zA(j)-105.37483D0)+1.4D0 !+1.5 Si es high tide
etaR=eta2(j)+0.18D0

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

zepA=(zt(Nx+3,j+2)-zt(Nx+2,j+2))/dep 
!zepA=(zt(Nx+2,j+2)-zt(Nx+1,j+2))/dep 

!--------------------------------
! Caracteristica Saliente
!Rmenos
h=q(1,:,j)	
u=q(2,:,j)
v=q(3,:,j)

!Solo tsunami
if (zA(j)>=zumbral) then !Borde Abierto
qA(1,j)=h(Nx)
qA(2,j)=u(Nx)
qA(3,j)=v(Nx)

else

if (h(Nx) == 0.0D0) then  !Caso seco
hL=0.0D0
uL=0.0D0
vL=0.0D0
RL=0.0D0


else
  !print*, j, borde
  call fzeroEP(j,borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epL,epxL,epyL,uL,vL,hL,zepL)
  !call fzeroEP_brent(borde,Fr2,Nx,dep,epx,epy,h,u,v,zep_x,dt,epL,epxL,epyL,uL,vL,hL,zepL)
  RL=(uL*epxL+vL*epyL)+2.0D0*sqrt(hL/Fr2*(epxL**2.0D0+epyL**2.0D0))
  !print*, RL, j
end if
     
if (fopt==0) then
tauL=0.0D0
Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
else
    if(hL/=0.0D0) then
    C=MC(Nx,j)
    call  tauU(tipo,C,Fr2,hL,uL,vL,tauL)
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))-tauL/hL*dt*(epxL+epyL)
    else
    tauL=0.0D0
    Rmas=RL-0.5D0*dt/Fr2*(zepL*(epxL**2.0D0+epyL**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if
end if


! Caracteristica Entrante

!----------------------------------
! En la perpendicular al borde
!qp=qx(j)*cos(alfa)+qy(j)*sin(alfa)

!----------------------------------
! Esto es solo para 8.6 y 8.8, por la subida de altura
!----------------------------------
! u0=qx(j)/(h0-1.4D0) !vuelvo al u original
! v0=qy(j)/(h0-1.4D0) !vuelvo al v original
! 
! qx2=u0*(h0+etaL)
! qy2=v0*(h0+etaL)

! En la direccion Xi
qp=qx(j)*epxA+qy(j)*epyA
qpar=qx(j)*epxA2+qy(j)*epyA2

! ! ---------------------------------------------
! Calculando RL con ambos datos desde afuera
!qpar=-qx(j)*sin(alfa)+qy(j)*cos(alfa)


 C0=sqrt(h0/Fr2)
epR=epA-dt*C0
hR=etaR+h0
! print*, 'L=',zA(j)-207.9032D0, 'h0=',h0, etaL, hL
! pause
if (hR<=0.0D0) then
hR=0.0D0
uR=0.0D0
vR=0.0D0
else

uR=-qx(j)!/hL  !SOLO para Columbia River , -0.5m/s
vR=qy(j)!/hL  !Solo para Columbia River

end if

call interp1(Nx,coordxi,zep_x,epR,zepR)
call interp1(Nx,coordxi,epx,epR,epxR)
call interp1(Nx,coordxi,epy,epR,epyR)

! Usando como perpendicular la direccion Xi
RR=(uR*epxR+vR*epyR)+2.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))

!Usando la recta perpendicular al borde
!RL=uL+2.0D0*sqrt(hL/Fr2)

!Rmas
if (fopt==0) then
tauR=0.0D0
else
  if (hR/=0.0D0) then
  C=MC(1,j)
  call tauU(tipo,C,Fr2,hR,uR,vR,tauR)
  tauR=tauR/hR
  else
  tauR=0.0D0
  end if
end if
Rmenos=RR-dt/Fr2*0.5D0*(zepR*epxR**2.0D0+zepA*epxA**2.0D0)-tauR*(epxR+epyR)*dt

!!!-Solucion CC---------

hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)

Uper=(Rmas+Rmenos)/2.0D0

UparL=uL*epxA2+vL*epyA2 !U2
UparR=uR*epxA2+vR*epyA2 !U2

if (Uper>0.0D0) then
Upar=UparL
! print*,'Upar=',Upar
else
Upar=UparR
! print*,'Upar=',Upar
end if

uA=Uper*epyA2/JacA+Upar*-epyA/JacA !U1*eta_y/J-Upar*Xiy/Jac
vA=Uper*-epxA2/JacA+Upar*epxA/JacA
! uA=Uper*epyA2/JacA+Upar*-epyA/JacA !U1*eta_y/J-Upar*Xiy/Jac
! vA=Uper*-epxA2/JacA+Upar*epxA/JacA


! end if


! ! ---------------------------------------
! ! Usando metodología de Brett Sanders
! ! ---------------------------------------
! --------------------
! ! ! Completo con hR, uR, vR el q0 en la primera iteracion
! ! call q0(borde,pasoRK,hL,uL,vL,j)
!-----------------
! ! alfa=angulo2(j)
! ! 
! ! !Solucion Ec. Orden 3
! ! hU=q0_global(1,Nx,j)
! ! uU=q0_global(2,Nx,j)
! ! vU=q0_global(3,Nx,j)
! ! 
! ! U1u=uU*epxA+vU*epyA
! ! RmasU=U1u+2.0D0*sqrt(hU/Fr2)*(epxA**2.0D0+epyA**2.0D0)
! ! 
! ! qp=qx(j)*cos(alfa)+qy(j)*sin(alfa)
! ! K=Fr2/(32.0D0*(epxA**2.0D0+epyA**2.0D0))
! ! a1=K
! ! a2=-K*RmasU
! ! a3=-K*RmasU**2.0D0
! ! a4=K*RmasU**3.0D0-qp
! ! call grado3(a1,a2,a3,a4,Rmenos)
! ! 
! ! hA=Fr2/16.0D0*(Rmas-Rmenos)**2.0D0/(epxA**2.0D0+epyA**2.0D0)
! ! U1=(Rmas+Rmenos)/2.0D0
! ! 
! ! uR=qx(j)/hA
! ! vR=qy(j)/hA
! ! UparL=-uL*sin(alfa)+vL*cos(alfa)
! ! UparR=-uR*sin(alfa)+vR*cos(alfa)
! ! ! UparL=uL*epxA2+vL*epyA2
! ! ! UparR=uR*epxA2+vR*epyA2
! ! 
! ! if (U1>0.0D0) then
! ! Upar=UparL
! ! else
! ! Upar=UparR
! ! end if
! ! 
! ! 
! ! !uA=U1*epyA2/JacA-Upar*epyA/JacA
! ! !vA=-U1*epxA2/JacA+Upar*epxA/JacA
! ! 
! ! uA=U1*cos(alfa)-Upar*sin(alfa)
! ! vA=U1*sin(alfa)+Upar*cos(alfa)

if (hA<1e-7) then
qA(1,j)=0.0D0
qA(2,j)=0.0D0
qA(3,j)=0.0D0
else
qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA
end if


end if !End if the zumbral

! print*, j, qA(:,j)!,uA, vA
! pause
END DO



END SUBROUTINE INFLOWN_xi
