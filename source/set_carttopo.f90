subroutine set_carttopo
  use mpi
  use mpi_surf
  use global_variables
  
  implicit none
  !split current zone in a cart topo
  !catch periodic boundaries here
  
  
  isperiodic(1)=.false.
  isperiodic(2)=.false.
  reorder=.true.  
  if (CB_real(1)==2.and.CB_real(2)==2) then
    isperiodic(1)=.true.
  else if (CB_real(1)==2.and.CB_real(2)/=2) then
    print*,'wrong boundary in bc(1) or bc(2) (periodic)'
    call mpi_abort(mpi_comm_world,100,ierror)
  else if (CB_real(1)/=2.and.CB_real(2)==2) then
    print*,'wrong boundary in bc(1) or bc(2) (periodic)'
    call mpi_abort(mpi_comm_world,100,ierror)
  end if
  
  if (CB_real(3)==2.and.CB_real(4)==2) then
    isperiodic(2)=.true.
  else if (CB_real(3)==2.and.CB_real(4)/=2) then
    print*,'wrong boundary in bc(3) or bc(4) (periodic)'
    call mpi_abort(mpi_comm_world,100,ierror)
  else if (CB_real(3)/=2.and.CB_real(4)==2) then
    print*,'wrong boundary in bc(3) or bc(4) (periodic)'
    call mpi_abort(mpi_comm_world,100,ierror)
  end if
  
  !mpi knows how to split comm_world better into a cart topo
  call mpi_dims_create(nproc,ndim,dims,ierror)
  call mpi_cart_create(mpi_comm_world,ndim,dims, &
	isperiodic,reorder,comm2d,ierror)
  call mpi_comm_rank(comm2d,myrank2d,ierror) 
  !get the coords for this core
  call mpi_cart_get(comm2d,ndim,dims,isperiodic,coords,ierror)
  
  !know my neighbors
  !left/right & back/front as when looking to the matrix
  !increasing indices to the upper right corner 
  !ex:  1,0	1,1
  !	0,0	0,1  
  !left(1,1)=1,0;  back(1,1)=0,1
  !right(0,0)=0,1; front(0,0)=1,0
  !left(1,0) =-2=mpi_proc_null
  !right(0,1)=-2=mpi_proc_null
  call mpi_cart_shift(comm2d, 0, shift, myback, myfront, ierror)
  call mpi_cart_shift(comm2d, 1, shift, myleft, myright, ierror)
   
end subroutine