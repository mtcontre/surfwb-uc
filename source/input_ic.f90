subroutine input_ic
  use mpi_surf
  use global_variables
  implicit none
  
  integer :: i,j
  allocate (q_buff(3,nxi,neta))
  
  select case (initqopt)
    case(1)
      open(unit=2, file=initqname(1))
      read(2,*) ((q_buff(1,i,j),j=1,neta),i=1,nxi)
      close(unit=2)
      
      open(unit=3, file=initqname(2))
      read(3,*) ((q_buff(2,i,j),j=1,neta),i=1,nxi)
      close(unit=3)
      
      open(unit=4, file=initqname(3))
      read(4,*) ((q_buff(3,i,j),j=1,neta),i=1,nxi)
      close(unit=4)  
  end select


end subroutine