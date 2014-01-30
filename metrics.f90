!Subrutina que crea las mallas xix, xiy, etax, etay
SUBROUTINE metrics 
  USE MPI
  USE geometries
  USE global_variables
  USE Jacobianos
  implicit none
  integer		::i,j,ierror
  real (kind=8)			::gj
  !real (kind=8)			::Sc,Se,Sz,SL


!Define tamaño de matrices
  allocate(xc(Nbx,Nby),yc(Nbx,Nby),&
  xe(Nbx,Nby),ye(Nbx,Nby)) 
  
  allocate(aj_global(Nbx,Nby),xi_global(2,Nbx,Nby),&
  eta_global(2,Nbx,Nby),Jac_global(Nbx,Nby),Jac_global_xi(Nbx+1,Nby),Jac_global_eta(Nbx,Nby+1))
  


!Primero calculo celdas centrales de 2 a N-1, xi_x y xi_y

!do nz=1,nzone
!Calculo celdas 2 a Nbx-1 para cada j
do j=1,Nby; do i=2,Nbx-1
	xc(i,j)=0.5D0*(x_global(i+1,j)-x_global(i-1,j))/dxi
	yc(i,j)=0.5D0*(y_global(i+1,j)-y_global(i-1,j))/dxi
end do; end do

!Calculo celdas 1 y Nbx para cada j 

do j=1,Nby
	xc(1,j)=(x_global(2,j)-x_global(1,j))/dxi
	yc(1,j)=(y_global(2,j)-y_global(1,j))/dxi
	xc(Nbx,j)=(x_global(Nbx,j)-x_global(Nbx-1,j))/dxi
	yc(Nbx,j)=(y_global(Nbx,j)-y_global(Nbx-1,j))/dxi
end do



do j=2,Nby-1; do i=1,Nbx
	xe(i,j)=0.5D0*(x_global(i,j+1)-x_global(i,j-1))/deta
	ye(i,j)=0.5D0*(y_global(i,j+1)-y_global(i,j-1))/deta
end do; end do


do i=1,Nbx
	xe(i,1)=(x_global(i,2)-x_global(i,1))/deta
	ye(i,1)=(y_global(i,2)-y_global(i,1))/deta
	xe(i,Nby)=(x_global(i,Nby)-x_global(i,Nby-1))/deta
	ye(i,Nby)=(y_global(i,Nby)-y_global(i,Nby-1))/deta
end do



! After calculating the values of xc,xe,yc,...etc, we
! calculate the metrics and the jacobian of the transformation:

do i=1,Nbx; do j=1,Nby
	gj=xc(i,j)*ye(i,j)-xe(i,j)*yc(i,j)
	if (abs(gj).le.1.0e-12) then
		print *,'zero jacobian in node=',i,j
		print *,gj
		print*, 'xc(i,j)=', xc(i,j)
		print*, 'ye(i,j)=', ye(i,j)
		print*, 'xe(i,j)=', xe(i,j)
		print*, 'yc(i,j)=', yc(i,j)
		call mpi_abort(mpi_comm_world,1000,ierror)
	end if
	aj_global(i,j)=1.0D0/gj !Jacobiano, es una varibale global, J(i,j)
	xi_global(1,i,j)=aj_global(i,j)*ye(i,j) !xi_x
	xi_global(2,i,j)=-aj_global(i,j)*xe(i,j) !xi_y
	eta_global(1,i,j)=-aj_global(i,j)*yc(i,j) !eta_x
	eta_global(2,i,j)=aj_global(i,j)*xc(i,j) !eta_y
end do; end do


call jacs(xc,xe,yc,ye,Nbx,Nby,Jac_global_xi,Jac_global_eta)

Jac_global=aj_global
! print*, eta_global(2,1,:)
! pause

END SUBROUTINE metrics


!Funcion que calcula las velocidades contravariantes tomando como input
!q=(h,u,v), xi_x, xi_y, eta_x, eta_y, im, jm que debe ser Nbx y Nby más los bordes
SUBROUTINE contra(qint,csi,eta,ucn,im,jm)
  implicit none
  integer	::im,jm
  real (kind=8),dimension(4,im,jm)	::qint
  real (kind=8),dimension(2,im,jm)	::ucn,csi,eta

  integer	::i,j
!print *,'q2 in contra',qint(2,1:3,15)
!print *,'q3 in contra',qint(3,1:3,15)

do j=1,jm; do i=1,im
        ucn(1,i,j)=csi(1,i,j)*qint(2,i,j)+csi(2,i,j)*qint(3,i,j)
        ucn(2,i,j)=eta(1,i,j)*qint(2,i,j)+eta(2,i,j)*qint(3,i,j)
end do; end do

END SUBROUTINE contra

SUBROUTINE jacs(xc,xe,yc,ye,im,jm,jac1,jac2)
  implicit none
  integer	::im,jm
  real (kind=8),dimension(im,jm)	::xc,yc,xe,ye
  real (kind=8),dimension(im+1,jm)	::jac1,G1,xcm1,ycm1,xem1,yem1
  real (kind=8),dimension(im,jm+1)	::jac2,G2,xcm2,ycm2,xem2,yem2
  real (kind=8),dimension(2,im,jm)	::xi,eta
  integer	::i,j

!Calculo Metricas intermedias
do i=1,im-1; do j=1,jm
xcm1(i+1,j)=(xc(i+1,j)+xc(i,j))/2.0D0
ycm1(i+1,j)=(yc(i+1,j)+yc(i,j))/2.0D0
xem1(i+1,j)=(xe(i+1,j)+xe(i,j))/2.0D0
yem1(i+1,j)=(ye(i+1,j)+ye(i,j))/2.0D0
end do; end do
do j=1,jm
xcm1(1,j)=xc(1,j)
xcm1(im+1,j)=xc(im,j)
ycm1(1,j)=yc(1,j)
ycm1(im+1,j)=yc(im,j)

xem1(1,j)=xe(1,j)
xem1(im+1,j)=xe(im,j)
yem1(1,j)=ye(1,j)
yem1(im+1,j)=ye(im,j)
end do

do i=1,im; do j=1,jm-1
xcm2(i,j+1)=(xc(i,j+1)+xc(i,j))/2.0D0
ycm2(i,j+1)=(yc(i,j+1)+yc(i,j))/2.0D0
xem2(i,j+1)=(xe(i,j+1)+xe(i,j))/2.0D0
yem2(i,j+1)=(ye(i,j+1)+ye(i,j))/2.0D0
end do; end do

do i=1,im
xcm2(i,1)=xc(i,1)
xcm2(i,jm+1)=xc(i,jm)
ycm2(i,1)=yc(i,1)
ycm2(i,jm+1)=yc(i,jm)

xem2(i,1)=xe(i,1)
xem2(i,jm+1)=xe(i,jm)
yem2(i,1)=ye(i,1)
yem2(i,jm+1)=ye(i,jm)
end do

!Jacobianos
do i=1,im+1; do j=1,jm
G1(i,j)=xcm1(i,j)*yem1(i,j)-xem1(i,j)*ycm1(i,j)
Jac1(i,j)=1.0D0/G1(i,j)
if (G1(i,j)==0.0D0) then
print*, i,j
stop
end if
end do; end do

do i=1,im; do j=1,jm+1
G2(i,j)=xcm2(i,j)*yem2(i,j)-xem2(i,j)*ycm2(i,j)
Jac2(i,j)=1.0D0/G2(i,j)
if (G2(i,j)==0.0D0) then
print*, i,j
stop
end if
end do; end do

END SUBROUTINE jacs
