subroutine init_seq
  use global_variables
  implicit none
  integer:: i,j
  real(kind=8) :: U1,U2
  call adimension
  
  call metrics
  
  ! para el cfl despues				!esto tb se deberia ordenar, por ej en: adimension-metrics-celerities
  allocate(S1_global(nbx,nby),S2_global(nbx,nby))
  do i=1,Nbx; do j=1,Nby
    U1=qold_global(2,i,j)*xi_global(1,i,j)+qold_global(3,i,j)*xi_global(2,i,j)
    U2=qold_global(2,i,j)*eta_global(1,i,j)+qold_global(3,i,j)*eta_global(2,i,j)
    S1_global(i,j)=abs(U1)+C_global(i,j)*sqrt(xi_global(1,i,j)**2+xi_global(2,i,j)**2)
    S2_global(i,j)=abs(U2)+C_global(i,j)*sqrt(eta_global(1,i,j)**2+eta_global(2,i,j)**2)
  end do; end do

end subroutine

