subroutine mkgrid(Nx,Ny,r0,r1,ds1,x1,y1)

! ifort -o mkgrid main.f90 /usr/local/tecplot/lib/tecio64.a -lstdc++

implicit none
integer i,j,k,m,im(6),jm(6),km(6), Nx, Ny
integer n,nzone,nz,i0, i1
real LX, LY, pi, r1, r0, alpha(300)
real ds1, ds2, dx, angle_x, r2
real radi, theta, psic, A

real, allocatable:: z(:)
real, allocatable:: x(:,:,:,:), y(:,:,:,:), x1(:,:), y1(:,:)
real, allocatable:: xtec(:,:), ytec(:,:),ztec(:,:)
real, allocatable:: xtop(:,:), ytop(:,:)
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
allocate( x(Nx,Ny,1,nzone), y(Nx,Ny,1,nzone) )


! O-grid,  NZ = 1
! - - - - - - - - - - - - - - - - - - - - 
im(1) = Nx	!207
jm(1) = Ny	!117
km(1) = 1

pi = 3.14159265451

!r0 = 0.50E+00   	   ! internal radius 
!r1 = 2.06E+00	   ! radius of the grid

ds1 = 2.0E+00*pi/real(im(1)-1)

! Internal
alpha(1)=0.0E+00
x(1,1,1,1)=r0*cos(alpha(1))
y(1,1,1,1)=-r0*sin(alpha(1))
do i=2,im(1)
	alpha(i)=alpha(i-1)+ds1
	x(i,1,1,1)=r0*cos(alpha(i))
	y(i,1,1,1)=-r0*sin(alpha(i))
end do


! External
alpha(1)=0.0E+00
x(1,jm(1),1,1)=r1*cos(alpha(1))
y(1,jm(1),1,1)=-r1*sin(alpha(1))
do i=2,im(1)
	alpha(i)=alpha(i-1)+ds1
	x(i,jm(1),1,1)=r1*cos(alpha(i))
	y(i,jm(1),1,1)=-r1*sin(alpha(i))
end do


! Filling in
!ds1 = 2.0E-04
alpha = 0.0E+00

do i=1,im(1)
	radi = sqrt( x(i,1,1,1)**2.0E+00 + y(i,jm(1),1,1)**2.0E+00 )
      print *,i,radi
	call tanhs(jm(1),radi,ds1,0.00E+00,alpha(1:jm(1)))

	do j=1,jm(1)-1
		x(i,j,1,1) = x(i,1,1,1) + ( x(i,jm(1),1,1)-x(i,1,1,1) )*alpha(j)/radi
		y(i,j,1,1) = y(i,1,1,1) + ( y(i,jm(1),1,1)-y(i,1,1,1) )*alpha(j)/radi
	end do
end do


!-------------------------------------------------------------------

print *,' '
print *,'END GRID'

!TECPLO FILE GRID!

	write(filename, fmt = '(a,i6.6)')'grid'		!En filename escribe grid: cambiar grid para otros archivos

    I = TecIni('Cylinder'//NULLCHR,			&   	! title of file
		'X, Y'//NULLCHR,	&   			! list of variables
		trim(filename)//'.plt'//NULLCHR,	&   	! output file name
		'.'//NULLCHR,						&
		Debug,								&
		VIsDouble)



do nz=1,nzone

allocate( xtec(im(nz),jm(nz)), ytec(im(nz),jm(nz)), ztec(im(nz),jm(nz)) )
xtec=0.0
ytec=0.0
ztec=0.0

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
		I   = TecDat(III,xtec,0)
		I   = TecDat(III,ytec,0)

deallocate(xtec,ytec,ztec)
end do !NZ
    I   = TecEnd()

close(1)

! open  (unit = 41, file = 'grid.dat', form='unformatted')
! do nz = 1,nzone
! 
!           write  (unit = 41) nz
! 
! 		  write  (unit = 41) im(nz),jm(nz),km(nz)
! 
!           write  (unit = 41) (((x(i,j,k,nz),i=1,im(nz)), &
!                                            j=1,jm(nz)),k=1,km(nz))
!           write  (unit = 41) (((y(i,j,k,nz),i=1,im(nz)), &
!                                            j=1,jm(nz)),k=1,km(nz))
! end do
! close(unit = 41)

!Writing into global variables
do i=1,Nx; do j=1,Ny
x1(i,j)=x(i,j,1,nz)
y1(i,j)=y(i,j,1,nz)
end do; end do


end subroutine mkgrid



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
real xd,ds1,ds2

real s(nm), ds01, ds02, a, del, be, ik, sk, us, snd

ik=float(nm-1)							! I=nm-1
if (ds2.ne.0) then
	! Stretching from both sides
	ds01=ds1/xd							! non-dimensional first space
	ds02=ds2/xd							! non-dimensional last space
	a=sqrt(ds02/ds01)					! A given by eq. (C6) of the appendix
	be=1/(ik*sqrt(ds01*ds02))			! B, I took I=nm-1 instead of nm-3			(?)
	call solve(del,be)							! finds delta
	
	do n=1,nm
		sk=n-1							! correspondent value of s, from 0 to I (from 0 to ik)
		! u(s) from eq. (C6)
		us=0.5*(1+(tanh(del*(sk/ik-0.5))/(tanh(del/2))))
		snd=us/(a+us*(1-a))				! non-dimensional psi
		s(n)=xd*snd						! Psi (output)
	end do

else
	! Stretching only from n=1, see eqs. (55) and (56) Thompson et al.
	ds01=ds1/xd
	be=1/(ik*ds01)
	call solve(del,be)

	do n=1,nm
		sk=n-1
		snd=1+(tanh(0.5*del*(sk/ik-1))/tanh(del/2))
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
real del,dx,be
!print *,'B= ',be
del=3.15 ! First approximation for delta
dx=15.00  ! Arbitrary value 
do while (abs(dx).ge.1.0E-05)
	!print *,del,dx
	f=(sinh(del)/del)-be					! Function f(del)=0
	df=(cosh(del)/del)-(sinh(del)/(del*del))	! Derivative dF/d(del)
	dx=f/df
	del=del-dx								! Newton step
end do

end subroutine



