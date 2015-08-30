subroutine locate(xp1, yp1, zp1, ii1, jj1, kk1, &
			u_vel1, v_vel1, w_vel1, &
            status1)

USE kdtree2_module
USE LagrangianPT
USE global_variables
USE geometries
implicit none

!~ integer :: nz,im,jm,km
integer :: status1

real (kind = 8):: xp1, yp1, zp1, u_vel1, v_vel1, w_vel1
integer :: ii1,jj1,kk1, lt

integer :: i,j,k, inside
real  (kind = 8), dimension(1:8,1:3) :: cell
real  (kind = 8), dimension(1:3) :: point
real  (kind = 8)  , dimension(1:6) :: dis

! allocatable local variables for kdtree search:
type (kdtree2), pointer    :: tree
type (kdtree2_result), allocatable :: results(:)
real  (kind = 8), dimension(:,:), allocatable :: xt

status1 = 1

iflag = 0	
	inside = 0

	point(1)=xp1; point(2)=yp1; point(3)=zp1
!~ print *, 'LOCATE 1 : point n_LPT', point, n_LPT
!~ print *, treeflag
if (treeflag.eq.0) then

	ibck = max(1,ibck);	ifwd = min(Nbx-1,ifwd)
	jbck = max(1,jbck);	jfwd = min(Nby-1,jfwd)
	kbck = 1;	kfwd = 1
	

	! if I have to search the entire zone
	! better use k-dimensional tree:
	!
	!if (ibck.eq.1.and.ifwd.eq.(im-1)) then

else
	allocate(results(2))
       	lt = Nbx*Nby*Nbz
       	allocate(xt(3,lt))

       	xt(1,1:lt) = reshape(x_global(1:Nbx,1:Nby), &
              			(/Nbx*Nby/))
       	xt(2,1:lt) = reshape(y_global(1:Nbx,1:Nby), &
              			(/Nbx*Nby/))
       	xt(3,1:lt) = reshape(z_LPT(1:Nbx,1:Nby,1), &
              			(/Nbx*Nby*Nbz/))
       	tree => kdtree2_create(xt, rearrange=.TRUE.,sort=.FALSE.)

	call kdtree2_n_nearest(tree, point, 1, results)

		k = int((results(1)%idx - 1)/ (Nby*Nbx)) + 1
                results(1)%idx = results(1)%idx - Nby*Nbx * (k - 1)
          	j = int((results(1)%idx - 1) / Nbx) + 1
          	results(1)%idx = results(1)%idx - Nbx * (j - 1)
          	i = results(1)%idx
print *,'KDTREE', i,j,k
		ibck = max(1,i-8);	ifwd = min(Nbx-1,i+8)
		jbck = max(1,j-8);	jfwd = min(Nby-1,j+8)
		kbck = max(1,k-8);	kfwd = min(Nbz-1,k+8)

		call kdtree2_destroy(tree)
		deallocate (xt,results)

end if
!~ print *,'LOCATE ibck jwd', ibck, ifwd, jbck,jfwd
	do j=jbck,jfwd; do i=ibck,ifwd
           cell(1,1) = x_global(i  ,j  )
           cell(1,2) = y_global(i  ,j  )           
           cell(2,1) = x_global(i+1,j  )
           cell(2,2) = y_global(i+1,j  )
           cell(3,1) = x_global(i+1,j+1)
           cell(3,2) = y_global(i+1,j+1)           
           cell(4,1) = x_global(i  ,j+1)
           cell(4,2) = y_global(i  ,j+1)           
           cell(5,1) = x_global(i  ,j  )
           cell(5,2) = y_global(i  ,j  )           
           cell(6,1) = x_global(i+1,j  )
           cell(6,2) = y_global(i+1,j  )           
           cell(7,1) = x_global(i+1,j+1)
           cell(7,2) = y_global(i+1,j+1)           
           cell(8,1) = x_global(i  ,j+1)
           cell(8,2) = y_global(i  ,j+1)
           
		   cell(1,3) = z_LPT(i  ,j  ,1)
		   cell(2,3) = z_LPT(i+1,j  ,1)
		   cell(3,3) = z_LPT(i+1,j+1,1)
		   cell(4,3) = z_LPT(i  ,j+1,1)
		   cell(5,3) = z_LPT(i  ,j  ,2)
		   cell(6,3) = z_LPT(i+1,j  ,2)
		   cell(7,3) = z_LPT(i+1,j+1,2)
		   cell(8,3) = z_LPT(i  ,j+1,2)

!~ print *, 'LOCATE cell, i, j', cell, i, j

	     call search_c(point, cell, inside, dis)

    if (inside.eq.1) then
			ii1 = i; jj1= j; kk1=k
!~ 			print *,'LOCATE : search_c, i j k', i,j,k
!~ 			print *,'LOCATE 1 : u_vel1 v_vel1 w_vel1', u_vel1, v_vel1, w_vel1
			call linear_2D(Nbx,Nby,Nbz,ii1,jj1,kk1,&
					u_vel1,v_vel1,w_vel1, &
					xp1,yp1,zp1)
			
!~ 			call linear(Nbx,Nby,Nbz,ii1,jj1,kk1,&
!~ 					u_vel1,v_vel1,w_vel1, &
!~ 					xp1,yp1,zp1)
			iflag = 1
!~ 			print *, 'LOCATE : point, inside, dis', point, inside, dis
!~ 			print *,'LOCATE : xp1 yp1 zp1', xp1, yp1, zp1
!~ 			print *,'LOCATE 2: u_vel1 v_vel1 w_vel1', u_vel1, v_vel1, w_vel1
			xpn(n_LPT) = xp1; ypn(n_LPT) = yp1 ; zpn(n_LPT) = zp1
			u_vel(n_LPT)=u_vel1; v_vel(n_LPT)=v_vel1; w_vel(n_LPT)=w_vel1
			k1(1)=u_vel1; k1(2)=v_vel1; k1(3)=w_vel1;
!~ 			print *,'LOCATE : k1', k1
			return
    end if ! inside

              
	
	end do; end do

end subroutine locate




! auxiliary subroutines
!
SUBROUTINE search_c(point,cell,inside,dis)
  implicit none
  real (kind = 8)::point(3),cell(8,3),dis(6)
  integer::inside,i,p1(6),p2(6),p3(6),p4(6)
  inside=0
!~   print *, 'SEARCH_C point', point
  ! i=const planes
  p1(1)=1
  p2(1)=5
  p3(1)=8
  p4(1)=4

  p1(2)=2
  p2(2)=3
  p3(2)=7
  p4(2)=6
  
  ! j=const planes
  p1(3)=1
  p2(3)=2
  p3(3)=6
  p4(3)=5

  p1(4)=4
  p2(4)=8
  p3(4)=7
  p4(4)=3
  
  ! k=const planes
  p1(5)=1
  p2(5)=4
  p3(5)=3
  p4(5)=2

  p1(6)=5
  p2(6)=6
  p3(6)=7
  p4(6)=8

  do i=1,6
!~   print *, 'SEARCH_C i', i
     ! find the distance of the point from the face defined by the four points given by cell
     call vector(cell(p1(i),:),cell(p2(i),:),cell(p3(i),:),cell(p4(i),:),  &
          & point,dis(i))          
     if(abs(dis(i))<1.0E-06) dis(i)=0.
     if(dis(i).lt.0.0)then
!~      if(dis(i).LE.0.0)then
        inside=0
        
        exit
     else
        inside=1
     end if

  end do
!~ print *, 'SEARCH_C dis , inside', dis, inside
END SUBROUTINE search_c

SUBROUTINE vector(p1,p2,p3,p4,point,dis)
! calculating the distace of the point from the plane of four points p1..p4
! by approximation using a=p3-p1 b=p4-p2 ab=a x b (cross product of a & b)
! c is the vector from the point til middle of the four points p1..p4 (intersection of vectors a & b
! c.(axb)/||axb|| gives the distance from the plane of the two vectors a & b
  implicit none
  real (kind=8)::dis,p1(3),p2(3),p3(3),p4(3),point(3)
  real (kind=8)::a1v,a2v,a3v,b1v,b2v,b3v,c1v,c2v,c3v,ab1,ab2,ab3

  a1v=p4(1)-p2(1)
  a2v=p4(2)-p2(2)
  a3v=p4(3)-p2(3)
  
  b1v=p3(1)-p1(1)
  b2v=p3(2)-p1(2)
  b3v=p3(3)-p1(3)

  c1v=point(1)-(p1(1)+p2(1)+p3(1)+p4(1))/4.
  c2v=point(2)-(p1(2)+p2(2)+p3(2)+p4(2))/4.
  c3v=point(3)-(p1(3)+p2(3)+p3(3)+p4(3))/4.

  ab1=a2v*b3v-a3v*b2v
  ab2=a3v*b1v-a1v*b3v
  ab3=a1v*b2v-a2v*b1v
!~ print *,'VECTOR ab cv', ab1, ab2, ab3, c1v, c2v, c3v, point
!~ print *,'VECTOR p(i)', p1, p2, p3,p4
  dis=(c1v*ab1+c2v*ab2+c3v*ab3)/sqrt(ab1**2+ab2**2+ab3**2)
END SUBROUTINE vector



subroutine linear(Nbx,Nby,Nbz,ii1,jj1,kk1,&
		u_vel1,v_vel1,w_vel1,&
        xp1, yp1, zp1)

USE geometries
USE LagrangianPT, ONLY : q_LPT, zero

implicit none
integer :: Nbx,Nby,Nbz
integer :: ii1,jj1,kk1

real (kind = 8)  :: u_vel1, v_vel1, w_vel1, d1
real (kind = 8)  :: xp1,yp1,zp1

real  (kind=8), dimension(1:8) :: px, py, pz, q1, q2, q3, con
real  (kind=8), dimension(1:8) :: qv1, qv2, qv3, qd1, qd2, qd3
real  (kind=8), dimension(1:8,1:8) :: matrix
integer :: i,j,k,n, idx, itask
integer, dimension(1:8) :: indx

idx=1
d1=0.0000
!~ print *, 'LINEAR ii jj kk', ii1, jj1, kk1
!~ print *, 'LINEAR xp1', xp1, yp1
do k=kk1,kk1+1; do j=jj1,jj1+1; do i=ii1,ii1+1

	px(idx) = x_global(i,j)
	py(idx) = y_global(i,j)  
	pz(idx) = z_LPT(i,j,1)
	q1(idx) = q_LPT(1,i,j)
	q2(idx) = q_LPT(2,i,j)
	q3(idx) = q_LPT(3,i,j)

    qd1(idx) = zero !Du(1,i,j)
	qd2(idx) = zero !Du(2,i,j)
	qd3(idx) = zero !Du(3,i,j)

    qv1(idx) = zero !vort(1,i,j)
	qv2(idx) = zero !vort(2,i,j)
	qv3(idx) = zero !vort(3,i,j)

	
	idx = idx + 1
end do; end do; end do
!~ print *, 'Q LINEAR', q1, q_LPT(:,i,j)
!~ print *, 'LINEAR q1', q1
do n=1,8
	matrix(n,1)=px(n)*py(n)*pz(n)
	matrix(n,2)=px(n)*py(n)
	matrix(n,3)=px(n)*pz(n)
	matrix(n,4)=py(n)*pz(n)
	matrix(n,5)=px(n)
	matrix(n,6)=py(n)
	matrix(n,7)=pz(n)
	matrix(n,8)=1.0000E+00
end do

call ludcmp(matrix,8,8,indx,d1)


do itask = 1,3
	if (itask.eq.1) then
		do n=1,8
			con(n) = q1(n)
		end do
	end if
	if (itask.eq.2) then
		do n=1,8
			con(n) = q2(n)
		end do
	end if
	if (itask.eq.3) then
		do n=1,8
			con(n) = q3(n)
		end do
	end if
!~ 	print *, 'LINEAR con', matrix, con
	call lubksb(matrix,8,8,indx,con)

	if (itask.eq.1) then
		u_vel1 = con(1)*xp1*yp1*zp1 + con(2)*xp1*yp1 + con(3)*xp1*zp1 +&
				con(4)*yp1*zp1 + con(5)*xp1 + con(6)*yp1 +&
				con(7)*zp1 + con(8)
!~ 				print *, 'LINEAR con', con, xp1*yp1
	end if
	if (itask.eq.2) then
		v_vel1 = con(1)*xp1*yp1*zp1 + con(2)*xp1*yp1 + con(3)*xp1*zp1 +&
				con(4)*yp1*zp1 + con(5)*xp1 + con(6)*yp1 +&
				con(7)*zp1 + con(8)
	end if
	if (itask.eq.3) then
		w_vel1 = con(1)*xp1*yp1*zp1 + con(2)*xp1*yp1 + con(3)*xp1*zp1 +&
				con(4)*yp1*zp1 + con(5)*xp1 + con(6)*yp1 +&
				con(7)*zp1 + con(8)
	end if

!~ print *, 'LINEAR VEOLICITY', u_vel1, v_vel1, w_vel1

end do !itask

!~ print *, 'LINEAR xp1 yp1 zp1', xp1, yp1, zp1
end subroutine linear


!-----------------------------------------------------------------------


subroutine linear_2D(Nbx,Nby,Nbz,ii1,jj1,kk1,&
		u_vel1,v_vel1,w_vel1,&
        xp1, yp1, zp1)

USE geometries
USE LagrangianPT, ONLY : q_LPT, zero

implicit none
integer :: Nbx,Nby,Nbz
integer :: ii1,jj1, kk1

real (kind = 8)  :: u_vel1, v_vel1, w_vel1
real (kind = 8)  :: xp1,yp1,zp1
real (kind = 8)  :: r1, s1, x_result, x00, x01,x10, x11

r1 = (xp1-x_global(ii1,jj1))/(x_global(ii1+1,jj1)-x_global(ii1,jj1))
s1 = (yp1-y_global(ii1,jj1))/(y_global(ii1,jj1+1)-y_global(ii1,jj1))
x_result=zero

call blend_102a (r1, s1, q_LPT(1,ii1,jj1), q_LPT(1,ii1+1,jj1), &
					q_LPT(1,ii1,jj1+1), q_LPT(1,ii1+1,jj1+1), x_result)
!~ print*, 'LINEAR_2D u', q_LPT(1,ii1,jj1), q_LPT(1,ii1+1,jj1), &
!~ 					q_LPT(1,ii1,jj1+1), q_LPT(1,ii1+1,jj1+1)
u_vel1=x_result

call blend_102a (r1, s1, q_LPT(2,ii1,jj1), q_LPT(2,ii1+1,jj1), &
					q_LPT(2,ii1,jj1+1), q_LPT(2,ii1+1,jj1+1), x_result)

v_vel1=x_result

call blend_102a (r1, s1, q_LPT(3,ii1,jj1), q_LPT(3,ii1+1,jj1), &
					q_LPT(3,ii1,jj1+1), q_LPT(3,ii1+1,jj1+1), x_result)

w_vel1=x_result


 end

! =====================================================
! -------------------------------------------------------------

	SUBROUTINE ludcmp(a,n,np,indx,d)
        INTEGER n,np,indx(n),NMAX
        real d,a(np,np),TINY
        PARAMETER (NMAX=8,TINY=1.0e-20)
        INTEGER i,imax,j,k
        real aamax,dum,sum,vv(NMAX)
        d=1.0000000
        do i=1,n
            aamax=0.00000000
            do j=1,n
              if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j))
            enddo
            if (aamax .eq. 0.00000000) print*,'singular matrix in ludcmp'
            vv(i)=1.0000000/aamax
        enddo
        do j=1,n
            do i=1,j-1
               sum=a(i,j)
               do k=1,i-1
                 sum=sum-a(i,k)*a(k,j)
               enddo
               a(i,j)=sum
            enddo
            aamax=0.000000000
            do i=j,n
               sum=a(i,j)
               do k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
               enddo
               a(i,j)=sum
               dum=vv(i)*abs(sum)
               if (dum .ge.aamax) then
                  imax=i
                  aamax=dum
               endif
             enddo
             if (j.ne.imax)then
                do k=1,n
                   dum=a(imax,k)
                   a(imax,k)=a(j,k)
                   a(j,k)=dum
                enddo
                d=-d
                vv(imax)=vv(j)
             endif
             indx(j)=imax
             if(a(j,j) .eq.0.0000000) a(j,j)=TINY
             if(j.ne.n)then
                dum=1.00000000/a(j,j)
                do i=j+1,n
                  a(i,j)=a(i,j)*dum
                enddo
               endif
             enddo
             return
             END
! -------------------------------------------------------------------

              SUBROUTINE lubksb(a,n,np,indx,b)
              INTEGER n,np,indx(n)
              real  a(np,np),b(n)
              INTEGER i,ii1,j,ll
              real sum
              ii1 = 0
              do i=1,n
                 ll=indx(i)
                 sum=b(ll)
                 b(ll)=b(i)
                 if (ii1.ne.0) then
                   do j=ii1,i-1
                       sum=sum-a(i,j)*b(j)
                   enddo
                 else if (sum.ne.0.00000) then
                   ii1=i
                 endif
                 b(i)=sum
               enddo
               do i=n,1,-1
                  sum=b(i)
                 do j=i+1,n
                   sum=sum-a(i,j)*b(j)
                 enddo
                 b(i)=sum/a(i,i)
               enddo
               return
               END 

subroutine blend_102a ( r, s, x00, x01, x10, x11, x )
!
!*******************************************************************************
!
!! BLEND_102 extends scalar point data into a square.
!
!
!  Diagram:
!
!    01------------11
!     |      .      |
!     |      .      |
!     |.....rs......|
!     |      .      |
!     |      .      |
!    00------------10
!
!  Formula:
!
!    Written in terms of R and S, the map has the form:
!
!      X(R,S) =
!               1     * ( + x00 )
!             + r     * ( - x00 + x10 )
!             + s     * ( - x00       + x01 )
!             + r * s * ( + x00 - x10 - x01 + x11 )
!
!    Written in terms of the coefficients, the map has the form:
!
!      X(R,S) = x00 * ( 1 - r - s + r * s )
!             + x01 * (         s - r * s )
!             + x10 * (     r     - r * s )
!             + x11 * (             r * s )
!
!             = x00 * ( 1 - r ) * ( 1 - s )
!             + x01 * ( 1 - r ) *       s
!             + x10 *       r   * ( 1 - s )
!             + x11 *       r           s
!
!    The nonlinear term ( r * s ) has an important role:
!
!      If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in
!      a plane, and the mapping is affine.  All the interpolated data 
!      will lie on the plane defined by the four corner values.  In 
!      particular, on any line through the square, data values at 
!      intermediate points will lie between the values at the endpoints.  
!
!      If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
!      not lie in a plane, and the interpolation map is nonlinear.  On
!      any line through the square, data values at intermediate points
!      may lie above or below the data values at the endpoints.  The
!      size of the coefficient of r * s will determine how severe this
!      effect is.
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates where an interpolated value is 
!    desired.  
!
!    Input, real X00, X01, X10, X11, the data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S).
!

  implicit none
!!
  real (kind=8) :: r, s, x, x00, x01, x10, x11
!
  x =             + x00 &
      + r *     ( - x00 + x10 ) & 
      + s *     ( - x00       + x01 ) &
      + r * s * ( + x00 - x10 - x01 + x11 )

  return
end
