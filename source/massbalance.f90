SUBROUTINE massbalance
  !Calcula el "volumen" adimensional de agua dentro del dominio computacional

  use global_variables
  use geometries
  use mpi_surf
  use mpi
  ! use multi
  implicit none

  integer	:: i,j
  real (kind=8):: masaT, volT, rho
  real (kind=8), allocatable, dimension(:,:)	:: masa, vol

  rho=1000.0D0
  masaT=0.0D0
  volT=0.0D0
  allocate(masa(Nbx,Nby),vol(Nbx,Nby))

  do i=1,Nbx; do j=1,Nby
		  
		  vol(i,j)=dxi*deta*qold_global(1,i,j)
		  masa(i,j)=vol(i,j)*rho
		  volT=volT+vol(i,j)
		  masaT=masaT+masa(i,j)
  end do; end do

  call mpi_allreduce(volT,volT,1,mpi_double_precision,mpi_sum,comm2d,ierror)
  if (treal==tinit) then  
    vol0=volT 
  end if
  
  !Porcentaje de volumen con respecto a la condicion inicial
  Pvol=volT/vol0*100.0D0

  !print *,'VolT mass balance', volT
  !pause
END SUBROUTINE massbalance
