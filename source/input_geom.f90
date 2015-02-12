subroutine input_geom
  use mpi_surf
  use global_variables
  implicit none
  
  integer :: i,j
  allocate (x_buff(nxi,neta), y_buff(nxi,neta), z_buff(nxi,neta))
  
  select case (batiopt)
    case(1)
      open(unit=2, file=batiname(1))
      read(2,*) ((x_buff(i,j),j=1,neta),i=1,nxi)
      close(unit=2)
      
      open(unit=3, file=batiname(2))
      read(3,*) ((y_buff(i,j),j=1,neta),i=1,nxi)
      close(unit=3)
      
      open(unit=4, file=batiname(3))
      read(4,*) ((z_buff(i,j),j=1,neta),i=1,nxi)
      close(unit=4)  
  end select


end subroutine