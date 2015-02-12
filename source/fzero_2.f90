!fzero function for 2nd and 4th RK step
!Fzero no está siendo usada, pero sirve para genabs solamente

SUBROUTINE fzero0_2(Fr2,N,dep,epx,h,u,zep,epxA,hA,uA,dt,epR,epxR,uR,hR,zepR)
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

real (kind=8):: Fr2, dt, dep, zepR,epR, ur, hr, epxR, f, fp, err, epR1, epR2, tol, epR0, epxR0, uR0, hR0, fo,&
		epxA,hA,uA
!real (kind=8),dimension(:),allocatable:: ep,epx,h,u,zep
real (kind=8),dimension(N):: ep, epx, h, u,zep
integer:: N, it, itmax, i

!allocate(ep(N),epx(N),h(N),u(N))

it=0
itmax=100000
tol=10e-8
err=1.0D0

ep(1)=0.5D0
DO i=2,N
ep(i)=ep(i-1)+dep
END DO
epR1=0.0D0
!Metodo Secante
epR0=1.0D0
call interp1(N,ep,epx,epR0,epxR0)
call interp1(N,ep,u,epR0,uR0)
call interp1(N,ep,h,epR0,hR0)

DO
IF (err<tol) EXIT

IF (it>itmax) THEN
print*,'Max it alcanzado en N-R en fzero 0_2'
print*,err,(err-tol),epR1
print*, f, fp
stop
END IF
!IF (it>itmax) EXIT

call interp1(N,ep,epx,epR1,epxR)
call interp1(N,ep,u,epR1,uR)
call interp1(N,ep,h,epR1,hR)

f=epR1-dt*0.5D0*((uR-sqrt(hR/Fr2))*epxR+(uA-sqrt(hA/Fr2))*epxA)
!fp=1.0D0!+dt*(uR-sqrt(hR/Fr2))*epxR !Esto tengo que arreglarlo!!
fo=epR0+dt*0.5D0*((uR0-sqrt(hR0/Fr2))*epxR0+(uA-sqrt(hA/Fr2))*epxA)
fp=(fo-f)/(epR0-epR1)	!Metodo Secante

epR2=epR1-f/fp
! print*, uR,hR,uA,hA
! print*, f, fo,epR1, epR0

err=abs(epR2-epR1)!/abs(epR1)

epR0=epR1 !Método Secante
epxR0=epxR
uR0=uR
hR0=hR

epR1=epR2
it=it+1

! pause
END DO

call interp1(N,ep,zep,epR1,zepR)
epR=epR1


END SUBROUTINE fzero0_2

SUBROUTINE fzeroN_2(Fr2,N,dep,epx,h,u,zep,epxA,hA,uA,dt,epL,epxL,uL,hL,zepL)
USE coords
implicit none

real (kind=8):: Fr2, dt, dep, epL, uL, hL, epxL, f, fp, err, epL1, epL2, tol,zepL, fo, epL0, uL0,hL0,epxL0, &
		epxA,hA,uA
!real (kind=8),dimension(:),allocatable:: ep,epx,h,u
real (kind=8),dimension(N):: ep,epx,h,u, zep
integer:: N, it, itmax, i

!allocate(ep(N),epx(N),h(N),u(N))

it=0
itmax=100000
tol=10e-7
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
!print*, epxL0,uL0,hL0
!call interp1(N,ep,zep,epL0,zepL0)


DO

IF (err<tol) EXIT

IF (it>itmax) THEN
print*,'Max it alcanzado en N-R en FzeroN2'
stop
END IF

call interp1(N,ep,epx,epL1,epxL)
call interp1(N,ep,u,epL1,uL)
call interp1(N,ep,h,epL1,hL)
!print*, ep(N),u(N),h(N),zep(N)
! print*, epxL, uL, hL
!  
! pause
if (hL<0.0D0) then
hL=0.0D0
uL=0.0D0
end if

f=epL1-N+dt*0.5D0*((uL+sqrt(hL/Fr2))*epxL+(uA+sqrt(hA/Fr2))*epxA)
!fp=1.0D0!+dt*(uL+sqrt(hL/Fr2))*epxL
fo=epL0-N+dt*0.5D0*((uL0+sqrt(hL0/Fr2))*epxL0+(uA+sqrt(hA/Fr2))*epxA)	!Método Secante
fp=(fo-f)/(epL0-epL1)	!Metodo Secante

epL2=epL1-(f/fp)
! print*, 'epL1=',epL1,'N=',N,'dt=',dt,'uL=',uL,'hL=',hL,'epxL',epxL
! print*, 'fo=',fo,'f=',f,'epL0=',epL0,'epL1=',epL1
! print*, 'fp=', fp
! print*, 'epL2=',epL2
! pause

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
END SUBROUTINE fzeroN_2