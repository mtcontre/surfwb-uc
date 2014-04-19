SUBROUTINE OUTFLOW0(fopt,tipo,C,Nx,Ny,Fr2,dep,hs,Ns,t,dt,q,zt,ep_x,qA,zA)
!Subrutina que fija una altura de agua en el borde 0 del dominio

!Crear una que fije la altura en el borde N del dominio

USE coords

!Fija una altura afuera del dominio en el borde 0, flujo unidireccional v=0.0

implicit none
real (kind=8):: C,t,dt,Fr2,epA, hR, uR,vR, epR, epxR,epyR, RR,dep, zepR, epxA,epyA,zepA,tauR,Rmenos, hb,hA,zmin,uA,vA, Rmas, U1, alpha
integer:: Nx,Ny,j,i, Ns, fopt, tipo
real (kind=8), dimension(3,Nx,Ny)::q
real (kind=8), dimension(Nx+4,Ny+4)::zt
real (kind=8), dimension(2,Nx,Ny)::ep_x
real (kind=8), dimension(3,Ny)::qA
real (kind=8), dimension(Ny)::zA
real (kind=8), dimension(Nx)::h,u,v,zep_x,epx,epy
real (kind=8), dimension(Ns,2)::hs


!Busco Altura fija, se entrega la cota de la superficie libre

IF ((t+dt)>maxval(hs(:,1))) THEN
!      hb=hs(Ns,2)
     print*, 'Faltan datos hidrograma'
     stop
ELSE
    call interp1(Ns,hs(:,1),hs(:,2),(t+dt),hb)
    !hb es la cota (hs+zmin)
END IF
epA=0.0D0
zmin=minval(zt(3,3:Ny+2))

DO j=1,Ny


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
if (h(1) == 0.0D0) then  !Caso seco
hR=0.0D0
uR=0.0D0
vR=0.0D0
RR=0.0D0
!Rmas=0.0D0
else
  call fzeroN(Fr2,Nx,dep,epx,h,u,zep_x,dt,epR,epxR,uR,hR,zepR)
  call interp1(Nx,coordxi,epy,epR,epyR)
  !print*, uL,hL,epxL, zepL
  RR=(uR*epxR+vR*epyR)-2.0D0*sqrt(hR/Fr2*(epxR**2.0D0+epyR**2.0D0))
end if
     
if (fopt==0) then
tauR=0.0D0
else
    if(hR/=0.0D0) then
    call  tauU(tipo,C,Fr2,hR,uR,vR,tauR)
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))-tauR/hR*dt*(epxR+epyR)
    else
    tauR=0.0D0
    Rmenos=RR-0.5D0*dt/Fr2*(zepR*(epxR**2.0D0+epyR**2.0D0)+zepA*(epxA**2.0D0+epyA**2.0D0))
    end if
end if


!Borde: Metricas y Zts y CI para ec. implicita
zA(j)=zt(3,j+2)

!Situacion Borde
!Informacion entrante: altura fija.

if ((hb).le.zA(j)) then !nodo seco
hA=0.0D0
uA=0.0D0
vA=0.0D0
else
hA=(hb)-zA(j)
Rmas=4.0D0*epxA*sqrt(hA/Fr2)+Rmenos
U1=(Rmas+Rmenos)/(2.0D0*epxA)
alfa=angulo1(j)
uA=U1*cos(alfa)
vA=U1*sin(alfa)

end if
qA(1,j)=hA
qA(2,j)=uA
qA(3,j)=vA

END SUBROUTINE OUTFLOW0

! 
! SUBROUTINE OUTFLOWN
! !Subrutina que fija una altura de agua en el borde 0 del dominio
! 
! !Crear una que fije la altura en el borde N del dominio
! 
! 
! END SUBROUTINE OUTFLOWN