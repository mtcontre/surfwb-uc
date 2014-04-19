program mkgrid
! ifort -o mkgrid main.f90 /usr/local/tecplot/lib/tecio64.a -lstdc++

implicit none
integer i,j,k,m,im(6),jm(6),km(6)
integer n,nzone,nz,i0, i1
real (kind = 8)  LX, LY, pi, r1, r0, alpha(300)
real (kind = 8) ds1, ds2, dx, angle_x, r2
real (kind = 8) radi, theta, psic, A
real (kind = 8), allocatable:: z(:)
real (kind = 8), allocatable:: x(:,:,:,:), y(:,:,:,:), x1(:,:), y1(:,:)
real (kind = 8), allocatable:: xtec(:,:), ytec(:,:),ztec(:,:)
real (kind = 8), allocatable:: xtop(:,:), ytop(:,:)
!-------------------------------------------------
character (len = 256) :: filename
! variables to enable writing of TecPlot binary (*.plt) files

  integer (kind = 4)           :: TecIni, TecDat, TecZne
  integer (kind = 4)           :: TecEnd
  integer (kind = 4)           :: VIsDouble = 0
  integer (kind = 4)           :: Debug = 1
  
  integer (kind = 4)           :: III

  character (len = 1)          :: nullchr = char(0)
!-----------------------------------------------------

nzone = 1
allocate( x(300,300,1,nzone), y(300,300,1,nzone), x1(300,300),y1(300,300) )


! O-grid,  NZ = 1
! - - - - - - - - - - - - - - - - - - - - 
im(1) = 51
jm(1) = 26
km(1) = 1

pi = 3.14159265451D0

r0 = 2.0D0!E+00   	   ! internal radius 
r1 = 25.0D0!E+00	   ! radius of the grid

ds1 = 2.0D0*pi/real(im(1)-1)

! Internal
alpha(1)=0.0D0!E+00
x(1,1,1,1)=r0*cos(alpha(1))
y(1,1,1,1)=-r0*sin(alpha(1))
do i=2,im(1)
	alpha(i)=alpha(i-1)+ds1
	x(i,1,1,1)=r0*cos(alpha(i))
	y(i,1,1,1)=-r0*sin(alpha(i))
end do


! External
alpha(1)=0.0D0
x(1,jm(1),1,1)=r1*cos(alpha(1))
y(1,jm(1),1,1)=-r1*sin(alpha(1))
do i=2,im(1)
	alpha(i)=alpha(i-1)+ds1
	x(i,jm(1),1,1)=r1*cos(alpha(i))
	y(i,jm(1),1,1)=-r1*sin(alpha(i))
	!print*, alpha(i)
	
end do
!pause

! Filling in
ds1 = 0.005D0!E-3
alpha = 0.0D0!E+00

do i=1,im(1)
	radi = sqrt( x(i,1,1,1)**2.0D0 + y(i,jm(1),1,1)**2.0D0 )
	
	!radi=sqrt((x(i,1,1,1)-x(i,jm(1),1,1))**2.0E+00+(y(i,1,1,1)-y(i,jm(1),1,1)**2.0E+00))
	print *,i,radi, x(i,1,1,1)
	
	call tanhs(jm(1),radi,ds1,0.0D0,alpha(1:jm(1)))
		
	
	do j=1,jm(1)-1
		x(i,j,1,1) = x(i,1,1,1) + ( x(i,jm(1),1,1)-x(i,1,1,1) )*alpha(j)/radi
		y(i,j,1,1) = y(i,1,1,1) + ( y(i,jm(1),1,1)-y(i,1,1,1) )*alpha(j)/radi
	end do
end do


!-------------------------------------------------------------------

print *,' '
print *,'END GRID'


	write(filename, fmt = '(a,i6.6)')'grid'

    I = TecIni('Cylinder'//NULLCHR,			&   ! title of file
		'X, Y'//NULLCHR,	&   ! list of variables
		trim(filename)//'.plt'//NULLCHR,	&   ! output file name
		'.'//NULLCHR,						&
		Debug,								&
		VIsDouble)



do nz=1,nzone

allocate( xtec(im(nz),jm(nz)), ytec(im(nz),jm(nz)), ztec(im(nz),jm(nz)) )
xtec=0.0D0
ytec=0.0D0
ztec=0.0D0

k=1
do i=1,im(nz); do j=1,jm(nz)
	xtec(i,j) = x(i,j,k,nz)
	ytec(i,j) = y(i,j,k,nz)
end do; end do

	I = TecZne('Zone'//NULLCHR,    &     
			im(nz),			           &
			jm(nz),			           &
			1, &
			'BLOCK'//NULLCHR,		   &
			NULLCHR//NULLCHR)
			! total number of points
		III = jm(nz) * im(nz)
			! write each variable
			! the last argument in the following calls indicates
			! the precision of the the variable,
			! 0 = single, 1 = double
		I   = TecDat(III,xtec,1)
		I   = TecDat(III,ytec,1)

deallocate(xtec,ytec,ztec)
end do !NZ
    I   = TecEnd()

close(1)

do i=1,im(nzone); do j=1,jm(nzone)
x1(i,j)=x(i,j,1,nzone)
y1(i,j)=y(i,j,1,nzone)
end do; end do

open  (unit = 41, file = 'gridX.dat', form='unformatted')
          write  (unit = 41) ((x1(i,j),i=1,im(nzone)),j=1,jm(nzone))
close(unit = 41)

open  (unit = 42, file = 'gridY.dat', form='unformatted')
          write  (unit = 42) ((y1(i,j),i=1,im(nzone)),j=1,jm(nzone))
close(unit = 42)

open  (unit = 43, file = 'gridZ.dat', form='unformatted')
          write  (unit = 43) ((0.0D0,i=1,im(nzone)),j=1,jm(nzone))
close(unit = 43)

print*, im(nzone), jm(nzone)

end program



subroutine tanhs (nm,xd,ds1,ds2,s)
! This subroutine introduces 'nm' points along a segment
! of total longitude 'xd', by using the hyperbolic tangent stretching 
! method (see Chap. VIII, Thompson et al. 1985).
! nm = number of points
! xd = total length of the segment
! ds1= predefined first space (n=1)	dimensional
! ds2= predefined last space (n=nm) dimensional, if ds2=0 the stretching is made only around n=1
! s  = Output s(n), corresponds to the non-dimensional arc length given by eq.(C5) of the appendix
!      multiplied by the total length, dimensional variable

! rewritten by C.E. June, 2003
integer nm
real (kind = 8) xd,ds1,ds2

real (kind = 8) s(nm), ds01, ds02, a, del, be, ik, sk, us, snd

ik=float(nm-1)							! I=nm-1
if (ds2.ne.0) then
	! Stretching from both sides
	ds01=ds1/xd							! non-dimensional first space
	ds02=ds2/xd							! non-dimensional last space
	a=sqrt(ds02/ds01)					! A given by eq. (C6) of the appendix
	be=1.0D0/(ik*sqrt(ds01*ds02))			! B, I took I=nm-1 instead of nm-3			(?)
	call solve(del,be)							! finds delta
	
	do n=1,nm
		sk=n-1.0D0							! correspondent value of s, from 0 to I (from 0 to ik)
		! u(s) from eq. (C6)
		us=0.50D0*(1.0D0+(tanh(del*(sk/ik-0.5D0))/(tanh(del/2.0D0))))
		snd=us/(a+us*(1.0D0-a))				! non-dimensional psi
		s(n)=xd*snd						! Psi (output)
	end do

else
	! Stretching only from n=1, see eqs. (55) and (56) Thompson et al.
	ds01=ds1/xd
	be=1/(ik*ds01)
	call solve(del,be)

	do n=1,nm
		sk=n-1.0D0
		snd=1.0D0+(tanh(0.5D0*del*(sk/ik-1.0D0))/tanh(del/2.0D0))
		s(n)=xd*snd						! Psi (output)
	end do

end if
end subroutine

subroutine solve(del,be)
! This subroutine finds the value of
! the stretching variable 'del' (delta)
! in the hyperbolic tangent method
! by Newton-Raphson's iterations
!( see eq. 51, Chap. VIII, Thompson et al. 1985)
! B=be is the known variable (eq. 50)

! rewritten by C.E.  June, 2003
!Double Precision del,dx,be
real (kind = 8) del,dx,be
!print *,'B= ',be
del=3.15D0 ! First approximation for delta
dx=15.00D0  ! Arbitrary value 
do while (abs(dx).ge.1.0E-05)
	!print *,del,dx
	f=(sinh(del)/del)-be					! Function f(del)=0
	df=(cosh(del)/del)-(sinh(del)/(del*del))	! Derivative dF/d(del)
	dx=f/df
	del=del-dx								! Newton step
end do

end subroutine



