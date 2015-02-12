!Fzero no está siendo usada, pero sirve para genabs solamente

SUBROUTINE fzero0(Fr2,N,dep,epx,h,u,zep,dt,epR,epxR,uR,hR,zepR)
USE coords
implicit none
!N=numero nodos
!dep= dxi o det ==1
!dt=dt paso RK (0.5dt oficial)
!h=alturas
!u=velocidades
!zep=pendiente fondo dz/dep
!epx= xix o etay
!ep=epsilon, vector de xi o etas


real (kind=8):: Fr2, dt, dep, zepR,epR, ur, hr, epxR, f, fp, err, epR1, epR2, tol, epR0, epxR0, uR0, hR0, fo
!real (kind=8),dimension(:),allocatable:: ep,epx,h,u,zep
real (kind=8),dimension(N):: ep, epx, h, u,zep
integer:: N, it, itmax, i

!allocate(ep(N),epx(N),h(N),u(N))

it=0
itmax=100000
tol=10e-5
err=1.0D0

ep(1)=0.5D0
DO i=2,N
ep(i)=ep(i-1)+dep
END DO
epR1=0.5D0
!Metodo Secante
epR0=1.0D0
call interp1(N,ep,epx,epR0,epxR0)
call interp1(N,ep,u,epR0,uR0)
call interp1(N,ep,h,epR0,hR0)

DO

IF (err<tol) EXIT

IF (it>itmax) THEN
print*,'Max it alcanzado en N-R en fzero 0'
print*,err,(err-tol),epR1
print*, f, fp
stop
END IF

!print*,'fzero0'
call interp1(N,ep,epx,epR1,epxR)
call interp1(N,ep,u,epR1,uR)
call interp1(N,ep,h,epR1,hR)

if (hR0<0.0D0) then
hR0=0.0D0
uR0=0.0D0
end if
if (hR<0.0D0) then
hR=0.0D0
uR=0.0D0
end if

! print*,epR1,epxR,uR,hR, h(1),u(1),fp,fo, hR0
! print*, 'holo'



f=epR1+dt*(uR-sqrt(hR/Fr2))*epxR
!fp=1.0D0!+dt*(uR-sqrt(hR/Fr2))*epxR !Esto tengo que arreglarlo!!
fo=epR0+dt*(uR0-sqrt(hR0/Fr2))*epxR0
fp=(fo-f)/(epR0-epR1)	!Metodo Secante
! print*, 'f',fp,fo, epR0, epR1
! pause

epR2=epR1-f/fp
err=abs(epR2-epR1)!/abs(epR1)

! print*,err, (err-tol)
! pause

epR0=epR1 !Método Secante
epxR0=epxR
uR0=uR
hR0=hR

epR1=epR2
it=it+1

! IF (err<tol) THEN
! print*, 'error menor tolerancia'
! END IF 

END DO

call interp1(N,ep,zep,epR1,zepR)
epR=epR1


END SUBROUTINE fzero0
!----------------------------------------------------------------------------

SUBROUTINE fzeroN(Fr2,N,dep,epx,h,u,zep,dt,epL,epxL,uL,hL,zepL)
USE coords
implicit none
real (kind=8):: Fr2, dt, dep, epL, uL, hL, epxL, f, fp, err, epL1, epL2, tol,zepL, fo, epL0, uL0,hL0,epxL0
!real (kind=8),dimension(:),allocatable:: ep,epx,h,u
real (kind=8),dimension(N):: ep,epx,h,u, zep
integer:: N, it, itmax, i

!allocate(ep(N),epx(N),h(N),u(N))

it=0
itmax=100000
tol=10e-5
err=1.0D0
ep(1)=0.5D0
DO i=2,N
ep(i)=ep(i-1)+dep
END DO
epL1=ep(N)-0.5D0
!Metodo Secante
epL0=ep(N)

call interp1(N,ep,epx,epL0,epxL0)
call interp1(N,ep,u,epL0,uL0)
call interp1(N,ep,h,epL0,hL0)
!call interp1(N,ep,zep,epL0,zepL0)
! print*, epL0, h(N),u(N)
! print*, epxL0,uL0,hL0

DO

IF (err<tol) EXIT

IF (it>itmax) THEN
print*,'Max it alcanzado en N-R en FzeroN'
stop
END IF

call interp1(N,ep,epx,epL1,epxL)
call interp1(N,ep,u,epL1,uL)
call interp1(N,ep,h,epL1,hL)
if (hL<0.0D0) then
hL=0.0D0
uL=0.0D0
end if

f=epL1-N+dt*(uL+sqrt(hL/Fr2))*epxL
fo=epL0-N+dt*(uL0+sqrt(hL0/Fr2))*epxL0	!Método Secante
fp=(fo-f)/(epL0-epL1)	!Metodo Secante

epL2=epL1-(f/fp)
err=abs(epL2-epL1)!/abs(epL1)

epL0=epL1 !Método Secante
epxL0=epxL
uL0=uL
hL0=hL

epL1=epL2

it=it+1


END DO

call interp1(N,ep,zep,epL1,zepL)
epL=epL1


END SUBROUTINE fzeroN

!----------------------------------------------------------------------------

SUBROUTINE fzeroH(Rmas,q,epx,hs0,zmin,zA,hA,uA)

!Funcion que calcula la altura de aguas en un punto en la seccion
!resuelve implicitamente: Rmas-(uA+2C)epx=0
USE global_variables
USE geometries
real (kind=8),dimension(Nby,2)::seccion
real (kind=8):: Rmas,q,h0,zmin,zA,h1,h2,err,f,fp,hA, epx, uA,uAo,Ao,A,Pm,Tsup,Rh,Dh,hs0,hs1,hs2,Zh
integer:: itf, itmaxf, N
real (kind=8),dimension(Nby)::hi
itf=0
itmaxf=100000
tol=10e-5

hs1=hs0*1.1D0
seccion(:,1)=y_global(Nbx,:)
seccion(:,2)=z_global(Nbx,:)

!hs es la altura de agua en el punto mas bajo de la seccion
!CB sirve para escalas iguales a 1.

if (q/=0.0D0) then

      DO

      IF (itf>itmaxf) THEN
      print*,'Max it alcanzado en N-R en FzeroH'
      END IF

      IF (itf>itmaxf) EXIT

      !hs0=h0+zA-zmin
      !hs1=h1+zA-zmin
      
	    
      call PropSecNat(seccion,Nby,hs0,Ao,Pm,Tsup,Rh,Dh)
      call PropSecNat(seccion,Nby,hs1,A,Pm,Tsup,Rh,Dh)
      uAo=-q/Ao !qb entra en el borde xi=Nx, u va hacia -x
      uA=-q/A
      h0=hs0+zmin-zA
      h1=hs1+zmin-zA
      
      fo=Rmas-(uAo+2.0D0*sqrt(h0/Fr2))*epx
      f=Rmas-(uA+2.0D0*sqrt(h1/Fr2))*epx

      fp=(fo-f)/(hs0-hs1)	!Metodo Secante
!       print*,Rmas,uAo,h0,uA,h1,epx
!       print*,fo,f,fp
!       pause
      hs2=hs1-(f/fp)

      err=abs(hs2-hs1)/abs(hs1)
      hs1=hs2
      itf=itf+1

      
      IF (err<tol) EXIT 

      END DO
      
      if (zA<(hs1+zmin)) then !nodo mojado
      hA=hs1+zmin-zA
      call PropSecNat(seccion,Nby,hs1,A,Pm,Tsup,Rh,Dh)
      uA=-q/A
      !print*, uA, hA, hs1,zmin,zA
      !pause
      else
      hA=0.0D0
      uA=0.0D0
      end if
else
      
      if (zA<(hs1+zmin)) then
      hA=hs1+zmin-zA
      uA=0.0D0
      else
      hA=0.0D0
      uA=0.0D0
      end if
!       print*, 'Holo'
!       pause
end if

END SUBROUTINE fzeroH
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

SUBROUTINE fzeroH_2(Rmas,q,epx,hs0,zmin,hA,uA)

!Funcion que calcula la altura de aguas en un punto en la seccion
!resuelve implicitamente: Rmas-(uA+2C)epx=0
USE global_variables
USE geometries
real (kind=8),dimension(Nby,2)::seccion
real (kind=8):: Rmas,q,h0,zmin,zA,h1,h2,err,f,fp,hA, epx, uA,uAo,Ao,A,Pm,Tsup,Tsupo,Rh,Dh,hs0,hs1,hs2,Zh
integer:: itf, itmaxf, N
real (kind=8),dimension(Nby)::hi
itf=0
itmaxf=100000
tol=10e-5

hs1=hs0*1.1D0
seccion(:,1)=y_global(Nbx,:)
seccion(:,2)=z_global(Nbx,:)

!hs es la altura de agua en el punto mas bajo de la seccion
!CB sirve para escalas iguales a 1.


      DO

      IF (itf>itmaxf) THEN
      print*,'Max it alcanzado en N-R en FzeroH'
      END IF

      IF (itf>itmaxf) EXIT

      !hs0=h0+zA-zmin
      !hs1=h1+zA-zmin
      
	    
      call PropSecNat(seccion,Nby,hs0,Ao,Pm,Tsupo,Rh,Dh)
      call PropSecNat(seccion,Nby,hs1,A,Pm,Tsup,Rh,Dh)
      uAo=q/Ao !qb entra en el borde xi=Nx, u va hacia -x
      uA=q/A
      h0=Ao/Tsupo
      h1=A/Tsup
      
      fo=Rmas-(uAo+2.0D0*sqrt(h0/Fr2))*epx
      f=Rmas-(uA+2.0D0*sqrt(h1/Fr2))*epx

      fp=(fo-f)/(hs0-hs1)	!Metodo Secante
!       print*,Rmas,uAo,h0,uA,h1,epx
!       print*,fo,f,fp
!       pause
      hs2=hs1-(f/fp)

      err=abs(hs2-hs1)/abs(hs1)
      hs1=hs2
      itf=itf+1

      
      IF (err<tol) EXIT 

      END DO
      
      hA=hs1
      call PropSecNat(seccion,Nby,hs1,A,Pm,Tsup,Rh,Dh)
      uA=q/A
      
END SUBROUTINE fzeroH_2
