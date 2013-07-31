!-----------------------------------------------------------------------
! Esta funcion se ocupa para los genabs solamente
! Hay otra en fzeroEP_1.f90 para el inflow y el outflow

SUBROUTINE fzero_B(borde,Fr2,N,dep,epx,h,u,zep,dt,epR,epxR,uR,hR,zepR) ! Brent method
!USE nrtype; USE nrutil, ONLY :nrerror
USE coords
IMPLICIT NONE
!N=numero nodos
!dep= dxi o det ==1
!dt=dt paso RK (0.5dt oficial)
!h=alturas
!u=velocidades
!zep=pendiente fondo dz/dep
!epx= xix o etay
!ep=epsilon, vector de xi o etas
REAL(kind=8) :: x1,x2
REAL(kind=8) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,xm
INTEGER :: N,itmax,iter,i,borde
REAL (kind=8):: Fr2,dt,dep,zepR,epR,ur,hr,epxR,f,fp,err1,epR1,epR2,tol,epR0,epxR0,uR0,hR0,fo
REAL (kind=8),dimension(N):: ep,epx,h,u,zep

!Using Brent’s method, find the root of a function func known to lie between x1 and x2.
!The root, returned as zbrent, will be refined until its accuracy is tol.
!Parameters: Maximum allowed number of iterations, and machine floating-point precision.
 iter=0
itmax=10000
tol=10e-5
err1=1.0D0

ep(1)=0.5D0
DO i=2,N
ep(i)=ep(i-1)+dep
END DO

if (borde==1.OR.borde==3) then !xi=0 o eta=0
epR1=0.5D0
epR0=1.0D0
else
!print*, 'holo'
epR1=ep(N)-0.5D0
epR0=ep(N)
end if

call interp1(N,ep,epx,epR0,epxR0)
call interp1(N,ep,u,epR0,uR0)
call interp1(N,ep,h,epR0,hR0)
call interp1(N,ep,epx,epR1,epxR)
call interp1(N,ep,u,epR1,uR)
call interp1(N,ep,h,epR1,hR)

if (borde==1.OR.borde==3) then 
fo=epR0+dt*(uR0-sqrt(hR0/Fr2))*epxR0
f=epR1+dt*(uR-sqrt(hR/Fr2))*epxR
else
fo=epR0-N+dt*(uR0+sqrt(hR0/Fr2))*epxR0
f=epR1-N+dt*(uR+sqrt(hR/Fr2))*epxR
end if

a=epR0
b=epR1
fa=fo
fb=f

IF ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) THEN
!print*,'f(a) and f(b) of the same sign'
END IF
c=b
fc=fb

do iter=1,itmax

if (hR0<0.0D0) then
hR0=0.0D0
uR0=0.0D0
end if
if (hR<0.0D0) then
hR=0.0D0
uR=0.0D0
end if

if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
!Rename a, b, c and adjust bounding interval d.
c=a 
fc=fa 
d=b-a
e=d
end if

if (abs(fc) < abs(fb)) then
a=b
b=c
c=a
fa=fb
fb=fc
fc=fa
end if

xm=0.5D0*(c-b)
if (abs(xm) <= tol .or. fb == 0.0) then
epR=b
EXIT
end if
if (abs(e) >= tol .and. abs(fa) > abs(fb)) then
s=fb/fa !Attempt inverse quadratic interpolation.
if (a == c) then
p=2.0D0*xm*s
q=1.0D0-s
else
q=fa/fc
r=fb/fc
p=s*(2.0D0*xm*q*(q-r)-(b-a)*(r-1.0D0))
q=(q-1.0D0)*(r-1.0D0)*(s-1.0D0)
end if
if (p > 0.0D0) q=-q !Check whether in bounds.
p=abs(p)
if (2.0D0*p < min(3.0D0*xm*q-abs(tol*q),abs(e*q))) then
e=d !Accept interpolation.
d=p/q
else
d=xm !Interpolation failed; use bisection.
e=d
end if
else !Bounds decreasing too slowly; use bisection.
d=xm 
e=d
end if
a=b !Move last best guess to a.
fa=fb
b=b+merge(d,sign(tol,xm), abs(d) > tol ) !Evaluate new trial root.

call interp1(N,ep,epx,b,epxR)
call interp1(N,ep,u,b,uR)
call interp1(N,ep,h,b,hR)

if (borde==1.OR.borde==3) then
f=b+dt*(uR-sqrt(hR/Fr2))*epxR
else
f=b-N+dt*(uR+sqrt(hR/Fr2))*epxR
end if


IF (iter==itmax) THEN
print*,'zbrent:exceeded maximum iterations'
END IF
end do


epR=b
END SUBROUTINE fzero_B