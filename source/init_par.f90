subroutine init_par
  use mpi
  use mpi_surf
  use global_variables
  implicit none
  
  !parametros 		??mejor inicializar en el modulo?
  g=9.812D0
  hmin=1.0D-7
  
  !inicializar variables	??mejor en init_params?
  it=0
  dt=0.0D0
  t=0.0D0
  treal=0.0D0
  
  call mpi_comm_rank(mpi_comm_world, myrank, ierror)
  call mpi_comm_size(mpi_comm_world,nproc,ierror)
  
  !read input parameters
  !temporarly, everyone reads everything
  !then each processor takes a piece of data
  !and deallocates the unneccesary
  
  call input_control 
   
  !read geometries
  call input_geom
  
  !intial condition
  call input_ic	
  
  !friction matrix
  call input_friction
  
  !initialize other parameteres from read data
  call init_params
  
  !decompose the domain: each proc gets a piece of data
  !according to a cartesian topology
  call set_carttopo
  
  call decomp_domain
  
  call decomp_bcs
  
  call saveout_parameters	!save grids parameters to files
  
end subroutine