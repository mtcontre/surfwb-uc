SUBROUTINE metrics8 
  USE geometries
  USE global_variables
  implicit none
  integer		::i,j
  real (kind=8)			::gj
  real (kind=8), allocatable, dimension (:,:):: radios, alphas
!real (kind=8)			::Sc,Se,Sz,SL

!Define tamaño de matrices
  allocate(xc(Nbx,Nby),yc(Nbx,Nby),&
  xe(Nbx,Nby),ye(Nbx,Nby)) 
  
  allocate(aj_global(Nbx,Nby),xi_global(2,Nbx,Nby),&
  eta_global(2,Nbx,Nby))
  
!Obtengo radio y angulos:
!xi=alpha, eta=radio

allocate (radios(Nbx,Nby), alphas(Nbx,Nby))

open	(unit=101, file ='radios.dat', form='unformatted')
read	(unit=101) ((radios(i,j),i=1,Nbx),j=1,Nby)
close(unit=101)

open	(unit=102, file ='alphas.dat', form='unformatted')
read	(unit=102) ((alphas(i,j),i=1,Nbx),j=1,Nby)
close(unit=102)

! ds1=2.0D0*pi/real(Nbx-1)
! alphai(1)=0.0D0
! 
! do i=2,Nbx
! alphai(i)=alpha(i-1)+ds1
! end do
! do j=1,Nby
! 
! radios(i,j)=sqrt(x_global(i,j)**2+y_global(i,j)**2)
! alphas(i,j)=
! print*, alphas(i,j)
! end do
! end do

! After calculating the values of xc,xe,yc,...etc, we
! calculate the metrics and the jacobian of the transformation:

do i=1,Nbx; do j=1,Nby

	gj=radios(i,j)
	xc(i,j)=-radios(i,j)*sin(alphas(i,j))
	xe(i,j)=cos(alphas(i,j))
	
	yc(i,j)=-radios(i,j)*cos(alphas(i,j))
	ye(i,j)=-sin(alphas(i,j))
	
	if (gj.le.1.0e-12) then
		print *,'zero jacobian in node=',i,j
		print *,gj
		print*, 'xc(i,j)=',xc(i,j)
		print*, 'ye(i,j)=', ye(i,j)
		print*, 'xe(i,j)=',xe(i,j)
		print*, 'yc(i,j)=',yc(i,j)
		stop
	end if
	gj=radios(i,j)
	aj_global(i,j)=1.0D0/gj
! 	xi_global(1,i,j)=aj_global(i,j)*ye(i,j) !xi_x
! 	xi_global(2,i,j)=-aj_global(i,j)*xe(i,j) !xi_y
! 	eta_global(1,i,j)=-aj_global(i,j)*yc(i,j) !eta_x
! 	eta_global(2,i,j)=aj_global(i,j)*xc(i,j) !eta_y
	xi_global(1,i,j)=1.0D0/radios(i,j)*ye(i,j)
	xi_global(2,i,j)=-1.0D0/radios(i,j)*xe(i,j) !xi_y
	eta_global(1,i,j)=-1.0D0/radios(i,j)*yc(i,j) !eta_x
	eta_global(2,i,j)=1.0D0/radios(i,j)*xc(i,j) !eta_y
	print*, radios(i,j), alphas(i,j)
	
end do; end do






END SUBROUTINE metrics8
