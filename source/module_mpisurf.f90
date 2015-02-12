module mpi_surf
  integer :: ierror
  integer :: myrank, nproc
  integer, parameter :: master = 0
  
  
  !input buffers
  real (kind=8), dimension(:,:), allocatable :: x_buff, y_buff, z_buff, fric_buff
  real (kind=8), dimension(:,:,:), allocatable :: q_buff
  integer, dimension(4) :: CB_real
  
  !cart topo
  logical, dimension(2) :: isperiodic
  logical :: reorder
  integer, parameter :: ndim = 2
  integer, dimension(2) :: dims, coords
  integer:: myrank2d, comm2d, myleft, myright, myback, myfront
  integer, parameter :: shift=1
  integer ::  si, ei, sj, ej
  
end module