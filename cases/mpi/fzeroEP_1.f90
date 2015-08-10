SUBROUTINE fzeroEP(j,borde,Fr2,N,dep,epx,epy,h,u,v,zep,dt,epL,epxL,epyL,uL,vL,hL,zepL)
USE coords
implicit none
real (kind=8):: Fr2, dt, dep, epL, uL,vL, hL, epxL, epyL, f, fp, err, epL1, epL2, tol,zepL, fo, epL0, uL0,hL0,vL0,epxL0,epyL0,U1,U0
!real (kind=8),dimension(:),allocatable:: ep,epx,h,u
real (kind=8),dimension(N):: ep,epx,epy,h,u,v,zep
integer:: N, it, itmax, i, borde,j

!allocate(ep(N),epx(N),h(N),u(N))

it=0
itmax=100000
tol=10e-5
err=1.0D0
ep(1)=0.5D0
DO i=2,N
ep(i)=ep(i-1)+dep
END DO

if (borde==1.OR.borde==3) then !xi=0 o eta=0
epL1=0.5D0
epL0=0.0D0
else
!print*, 'holo' !Esto podria ser al reves, epL1=ep(N) epL0=ep(N)-0.5
epL1=ep(N)+0.5D0
epL0=ep(N)
end if


call interp1(N,ep,epx,epL0,epxL0)
call interp1(N,ep,epy,epL0,epyL0)
call interp1(N,ep,u,epL0,uL0)
call interp1(N,ep,v,epL0,vL0)
call interp1(N,ep,h,epL0,hL0)

! print*,'borde=',borde
! IF (borde==2) then
! 
! print*, 'hL0=',hL0
! print*, 'uL0=',uL0
! print*, 'vL0=',vL0
! print*, 'epxL0=',epxL0
! print*, 'epyL0=',epyL0
! 
! end if


!Metodo Secante
DO

IF (err<tol) EXIT

IF (it>itmax) THEN
print*,'Max it alcanzado en N-R en FzeroEP'
print*, borde, j
print*, 'holo'
print*, U1,uL,hL, epxL, vL, epyL
print*, U0,uL0,hL0, epxL0, vL0, epyL0
print*,'epL=',epL0,epL1,epL2
stop
END IF

call interp1(N,ep,epx,epL1,epxL)
call interp1(N,ep,epy,epL1,epyL)
call interp1(N,ep,u,epL1,uL)
call interp1(N,ep,v,epL1,vL)
call interp1(N,ep,h,epL1,hL)

! if (borde==2) then
! print*, 'hL=',hL
! print*, 'uL=',uL
! print*, 'vL=',vL
! print*, 'epxL=',epxL
! print*, 'epyL=',epyL
! print*, 'err=', err
! pause
! end if


if (hL0<0.0D0) then
hL0=0.0D0
uL0=0.0D0
vL0=0.0D0
end if


if (hL<0.0D0) then
hL=0.0D0
uL=0.0D0
vL=0.0D0
end if

U1=uL*epxL+vL*epyL
U0=uL0*epxL0+vL0*epyL0
!print*, 'U0=',U0

if (borde==1.OR.borde==3) then !xi=0 o eta=0
!print*, borde
f=epL1+dt*(U1-sqrt(hL/Fr2))*sqrt(epxL**2.0D0+epyL**2.0D0)
fo=epL0+dt*(U0-sqrt(hL0/Fr2))*sqrt(epxL0**2.0D0+epyL0**2.0D0)

else 
f=epL1-N+dt*(U1+sqrt(hL/Fr2))*sqrt(epxL**2.0D0+epyL**2.0D0)
fo=epL0-N+dt*(U0+sqrt(hL0/Fr2))*sqrt(epxL0**2.0D0+epyL0**2.0D0)

end if

!Metodo Secante

fp=(fo-f)/(epL0-epL1)	

epL2=epL1-(f/fp)

err=abs(epL2-epL1)!/abs(epL1)


! print*, f, fo
! print*, uL,U1,hL,epxL,epL1
! print*, uL0,U0,hL0,epxL0, epL0,epL2
! pause


epL0=epL1 !Método Secante
epxL0=epxL
uL0=uL
vL0=vL
hL0=hL

epL1=epL2

it=it+1


END DO

call interp1(N,ep,zep,epL1,zepL)
epL=epL1


END SUBROUTINE fzeroEP

!-----------------------------------------------------------------------------
!Fzero Metodo de Brent

SUBROUTINE fzeroEP_brent(borde,Fr2,N,dep,epx,epy,h,u,v,zep,dt,epL,epxL,epyL,uL,vL,hL,zepL)
USE coords
implicit none
real (kind=8) :: x1,x2
real (kind=8) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,xm
real (kind=8):: Fr2, dt, dep, epL, uL,vL, hL, epxL, epyL, f, fp, err1, epL1, epL2, tol,zepL, fo, epL0, uL0,hL0,vL0,epxL0,epyL0,U1,U0
!real (kind=8),dimension(:),allocatable:: ep,epx,h,u
real (kind=8),dimension(N):: ep,epx,epy,h,u,v,zep
integer:: N, iter, itmax, i, borde

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
epL1=0.5D0
epL2=0.0D0
else
epL1=ep(N)
epL0=ep(N)+0.5D0
end if

call interp1(N,ep,epx,epL0,epxL0)
call interp1(N,ep,epy,epL0,epyL0)
call interp1(N,ep,u,epL0,uL0)
call interp1(N,ep,v,epL0,vL0)
call interp1(N,ep,h,epL0,hL0)

call interp1(N,ep,epx,epL1,epxL)
call interp1(N,ep,epy,epL1,epyL)
call interp1(N,ep,u,epL1,uL)
call interp1(N,ep,v,epL1,vL)
call interp1(N,ep,h,epL1,hL)

U1=uL*epxL+vL*epyL
U0=uL0*epxL0+vL0*epyL0

!print*, 'U0=',U0

if (borde==1.OR.borde==3) then !xi=0 o eta=0
!print*, borde
f=epL1+dt*(U1-sqrt(hL/Fr2))*sqrt(epxL**2.0D0+epyL**2.0D0)
fo=epL0+dt*(U0-sqrt(hL0/Fr2))*sqrt(epxL0**2.0D0+epyL0**2.0D0)

else 
f=epL1-N+dt*(U1+sqrt(hL/Fr2))*sqrt(epxL**2.0D0+epyL**2.0D0)
fo=epL0-N+dt*(U0+sqrt(hL0/Fr2))*sqrt(epxL0**2.0D0+epyL0**2.0D0)

end if

a=epL0
b=epL1
fa=fo
fb=f

IF ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) THEN
!print*,'f(a) and f(b) of the same sign'
END IF
c=b
fc=fb

do iter=1,itmax

if (hL0<0.0D0) then
hL0=0.0D0
uL0=0.0D0
vL0=0.0D0
end if
if (hL<0.0D0) then
hL=0.0D0
uL=0.0D0
vL=0.0D0
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
epL=b
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

call interp1(N,ep,epx,b,epxL)
call interp1(N,ep,epy,b,epyL)
call interp1(N,ep,u,b,uL)
call interp1(N,ep,v,b,vL)
call interp1(N,ep,h,b,hL)

!f=b+dt*(uR-sqrt(hR/Fr2))*epxR

if (borde==1.OR.borde==3) then !xi=0 o eta=0
f=b+dt*(U1-sqrt(hL/Fr2))*sqrt(epxL**2.0D0+epyL**2.0D0)
else 
f=b-N+dt*(U1+sqrt(hL/Fr2))*sqrt(epxL**2.0D0+epyL**2.0D0)
end if


IF (iter==itmax) THEN
print*,'zbrent:exceeded maximum iterations'
END IF
end do

epL=b
call interp1(N,ep,zep,epL,zepL)

END SUBROUTINE fzeroEP_brent


SUBROUTINE fzeroEP_per(borde,Fr2,N,dep,epx,epy,h,u,v,zep,dt,epL,epxL,epyL,uL,vL,hL,zepL)
USE coords
implicit none
real (kind=8):: Fr2, dt, dep, epL, uL,vL, hL, epxL, epyL, f, fp, err, epL1, epL2, tol,zepL, fo, epL0, uL0,hL0,vL0,epxL0,epyL0,U1,U0
!real (kind=8),dimension(:),allocatable:: ep,epx,h,u
real (kind=8),dimension(N):: ep,epx,epy,h,u,v,zep
integer:: N, it, itmax, i, borde

!allocate(ep(N),epx(N),h(N),u(N))

it=0
itmax=100000
tol=10e-5
err=1.0D0
ep(1)=0.5D0
DO i=2,N
ep(i)=ep(i-1)+dep
END DO

if (borde==1.OR.borde==3) then !xi=0 o eta=0
epL1=0.5D0
epL0=0.0D0
else
!print*, 'holo' !Esto podria ser al reves, epL1=ep(N) epL0=ep(N)-0.5
epL1=ep(N)+0.5D0
epL0=ep(N)
end if


call interp1(N,ep,epx,epL0,epxL0)
call interp1(N,ep,epy,epL0,epyL0)
call interp1(N,ep,u,epL0,uL0)
call interp1(N,ep,v,epL0,vL0)
call interp1(N,ep,h,epL0,hL0)


!Metodo Secante
DO

IF (err<tol) EXIT

IF (it>itmax) THEN
print*,'Max it alcanzado en N-R en FzeroEP'
print*, borde
print*, 'holo'
print*, U1,uL,hL, epxL, vL, epyL
print*, U0,uL0,hL0, epxL0, vL0, epyL0
print*,'epL=',epL0,epL1,epL2
stop
END IF

call interp1(N,ep,epx,epL1,epxL)
call interp1(N,ep,epy,epL1,epyL)
call interp1(N,ep,u,epL1,uL)
call interp1(N,ep,v,epL1,vL)
call interp1(N,ep,h,epL1,hL)

if (hL0<0.0D0) then
hL0=0.0D0
uL0=0.0D0
vL0=0.0D0
end if


if (hL<0.0D0) then
hL=0.0D0
uL=0.0D0
vL=0.0D0
end if

! U1=uL*epxL+vL*epyL
! U0=uL0*epxL0+vL0*epyL0
U1=uL
U0=uL0
!print*, 'U0=',U0

if (borde==1.OR.borde==3) then !xi=0 o eta=0
!print*, borde
f=epL1+dt*(U1-sqrt(hL/Fr2))!*sqrt(epxL**2.0D0+epyL**2.0D0)
fo=epL0+dt*(U0-sqrt(hL0/Fr2))!*sqrt(epxL0**2.0D0+epyL0**2.0D0)

else 
f=epL1-N+dt*(U1+sqrt(hL/Fr2))!*sqrt(epxL**2.0D0+epyL**2.0D0)
fo=epL0-N+dt*(U0+sqrt(hL0/Fr2))!*sqrt(epxL0**2.0D0+epyL0**2.0D0)

end if

!Metodo Secante

fp=(fo-f)/(epL0-epL1)	

epL2=epL1-(f/fp)

err=abs(epL2-epL1)!/abs(epL1)

! print*, 'holi'
! print*, f, fo
! print*, uL,U1,hL,epxL,epL1
! print*, uL0,U0,hL0,epxL0, epL0,epL2
! pause


epL0=epL1 !Método Secante
epxL0=epxL
uL0=uL
vL0=vL
hL0=hL

epL1=epL2

it=it+1


END DO

call interp1(N,ep,zep,epL1,zepL)
epL=epL1


END SUBROUTINE fzeroEP_per






























