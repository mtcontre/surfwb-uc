subroutine decomp2d
  USE global_variables
  USE geometries
  USE MPI
  USE MPI_SURF
  

  !mpi knows how to split comm_world better into a cart topo:  
  !obtain dims vector given ndim and nproc (let mpi decide it)
  call mpi_dims_create(nproc,ndim,dims,ierror)
  !create the communicator comm2d
  call mpi_cart_create(mpi_comm_world,ndim,dims,isperiodic,&
    reorder,comm2d,ierror)
  !obtain myrank2d
  call mpi_comm_rank(comm2d,myrank2d,ierror)
  !get the topology_coords in the cart topo for this process
  call mpi_cart_get(comm2d,ndim,dims,isperiodic,topology_coords,ierror)  
  
  print*,'world rank =',myrank, '2d rank=', myrank2d, ' cartcoords=', &
    topology_coords(1),topology_coords(2)
  stop
!   call mpi_finalize(ierror)
end subroutine