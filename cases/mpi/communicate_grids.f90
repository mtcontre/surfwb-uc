subroutine communicate_grids(qbuff_recv,qbuff_send,bc,N)
  use mpi
  use multigrid_surf
  use global_variables
  implicit none
  integer :: N !length of this boundary
  integer :: left,right,back,front,shift
  integer ::bc !number of bc asked, if xi1,xiN,eta1 or etaN
  real(kind=8),dimension(N,2),allocatable::qbuff_recv,qbuff_send !buffers for send_recv operations
  
  shift=1
  !get my neighbour
  call mpi_cart_shift(comm2d,0,shift,back,front,ierror)
  
  call mpi_sendrecv(qbuff_send,2*N,mpi_double_precision,front,)
  MPI_SENDRECV(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype,
source, recvtag, comm, status)
call mpi_cart_shift(comm2d,1,shift,left,right,ierror)
end subroutine communicate_grids