SUBROUTINE ADIMENSION
  !Function that aplies the adimensionalization to the initial conditions
  !Funcion que aplica la adimensionalizacion a las condiciones inciales y a todo

  USE global_variables
  USE geometries
  implicit none
  integer :: i,j
  real (kind=8)::U1,U2
  allocate(V_global(nbx,nby), C_global(nbx,nby), VC(nbx,nby))
  do i=1,Nbx; do j=1,Nby
    x_global(i,j)=x_global(i,j)/L
    y_global(i,j)=y_global(i,j)/L
    z_global(i,j)=z_global(i,j)/H
    qold_global(1,i,j)=qold_global(1,i,j)/H
    qold_global(2,i,j)=qold_global(2,i,j)/U
    qold_global(3,i,j)=qold_global(3,i,j)/U			
    V_global(i,j)=sqrt((qold_global(2,i,j))**2.0D0+(qold_global(3,i,j))**2.0D0)			
    C_global(i,j)=sqrt(qold_global(1,i,j)/FR2)
    VC(i,j)=V_global(i,j)+C_global(i,j)			
  end do; end do
END SUBROUTINE ADIMENSION
