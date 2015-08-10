SUBROUTINE LPT ()

USE global_variables
USE LagrangianPT
USE geometries

implicit none

integer	:: i, j
!~ real (kind=8)	:: ind, dec, time
!~ logical	:: lexist, lexistT, lexistP, lexistT_eta
!~ character(len=1100):: filename, filenameT,filenameP,filenameT_eta, path
!~ character(len=10)::number, numbercaso
!~ character(len=1000)::intchar, ncaso
!~ real (kind=8), dimension(4)::x0, y0, i0, j0, H0
!~ real (kind=8), dimension(1)::m1,m2
!~ real (kind=8), dimension(4)::r, s, x, x00, x01, x10, x11 ! Interpolation by blending (see blend.f90, function blend_102)
!~ x0=(/5.01, 10.31, 15.315, 20.345 /)
!~ y0=(/8.1, 8.1, 8.1, 8.1 /)

num=500 ! number of particle
!~ print *, 'LPT it', it

IF (it.EQ.0.0D0) THEN ! Initializacion

allocate(xp(num), yp(num), zp(num), up(num),vp (num), wp(num))
allocate(xpn(num),ypn(num),zpn(num),upn(num),vpn(num),wpn(num))
allocate(u_vel(num),v_vel (num), w_vel(num))
allocate(ii(num),jj(num),kk(num),status_n(num))
allocate(z_LPT(Nbx,Nby,2))
allocate(q_LPT(3,Nbx,Nby))

!~ print *, 'xp iter', xpn, it

!DO i=1,num
DO i=1,25; DO j=1,20;

!	xp(i)=16.0D0
!	yp(i)=5.0D0+0.01D0*i

	xp(i+(j-1)*25)=19.0D0+(j-1)*0.1D0
	yp(i+(j-1)*25)=9.0D0 + i*0.08D0

!!~ 	xp(i+(j-1)*10)=15.0D0+(j-1)*0.5D0
!!~ 	yp(i+(j-1)*10)=7.0D0 + 0.5D0*i
	
	ii(i)=minloc(x_global(:,1)-xp(i),1,mask=(x_global(:,1)-xp(i)).ge.0)
	jj(i)=minloc(y_global(1,:)-yp(i),1,mask=(y_global(1,:)-yp(i)).ge.0)
	kk(i)=1

END DO; END DO
	

!DO i=1,250; DO j=1,2;
!
!!	xp(i)=16.0D0
!!	yp(i)=5.0D0+0.01D0*i
!
!	xp(i+(j-1)*250)=15.0D0+(j-1)*5.0D0
!	yp(i+(j-1)*250)=3.0D0 + i*0.1D0
!
!!!~ 	xp(i+(j-1)*10)=15.0D0+(j-1)*0.5D0
!!!~ 	yp(i+(j-1)*10)=7.0D0 + 0.5D0*i
!	
!	ii(i)=minloc(x_global(:,1)-xp(i),1,mask=(x_global(:,1)-xp(i)).ge.0)
!	jj(i)=minloc(y_global(1,:)-yp(i),1,mask=(y_global(1,:)-yp(i)).ge.0)
!	kk(i)=1
!
!
!END DO; END DO
	
zp=zero
up=zero
vp=zero
wp=zero
xpn=zero
ypn=zero
zpn=zero
upn=zero
wpn=zero
vpn=zero
status_n=1
kk=1.0D0

u_vel=zero
v_vel=zero
w_vel=zero
ibck=0
ifwd=0
jbck=0
jfwd=0
kbck=0
kfwd=0
k1=zero

z_LPT(:,:,1)=zero
z_LPT(:,:,2)=one

LPT_init=1

!~ print *, 'LPT ii jj', ii, jj

k1(3:6)=zero
END IF
! print *, 'LPT xp yp ', xp, yp
!pause
IF (LPT_init.EQ.1) THEN

LPT_init=0 ! first step
print *, 'LPT init OK'
ELSE

	call lpt_p()
	call lpt_c()

END IF
!~ print *, 'k1', k1

END SUBROUTINE
