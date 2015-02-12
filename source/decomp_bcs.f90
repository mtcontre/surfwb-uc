subroutine decomp_bcs
  use mpi
  use mpi_surf
  use global_variables
  
  implicit none
  
  CB = -1

  if (myback==mpi_proc_null) then
    CB(1) = CB_real(1)
  end if
  
  if (myfront==mpi_proc_null) then
    CB(2) = CB_real(2)
  end if
  
  if (myleft==mpi_proc_null) then
    CB(3) = CB_real(3)
  end if
  
  if (myright==mpi_proc_null) then
    CB(4) = CB_real(4)
  end if
  
  !periodic boundaries
  if (CB(1)==2) then
    CB(1)=-1
    CB(2)=-1
  end if
  
  if (CB(3)==2) then
    CB(3)=-1
    CB(4)=-1    
  end if

end subroutine