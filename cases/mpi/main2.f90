program main2
  use mpi
  use mpi_surf
  use global_variables

  implicit none
  
  call mpi_init(ierror)
    call init
    
  call mpi_finalize(ierror)
end program main2