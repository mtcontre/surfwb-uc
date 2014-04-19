subroutine bcast_bcs
  !set CB vector and scatter input time series if neccesary
  use global_variables
  use couplingbc
  use mpi_surf
  use multigrid_surf
  use mpi
  implicit none
  integer i
  integer, dimension(mpi_status_size) :: my_status
  integer,dimension(:),allocatable :: isend_request
  integer myrank_temp
  
  CB(1)=-1
  CB(2)=-1
  CB(3)=-1
  CB(4)=-1  
  
  if (myback==mpi_proc_null) then 
    CB(1)=CB_real(1)
  end if
  
  if (myfront==mpi_proc_null) then
    CB(2)=CB_real(2)
  end if
  
  if (myleft==mpi_proc_null) then
    CB(3)=CB_real(3)
  end if
  
  if (myright==mpi_proc_null) then
    CB(4)=CB_real(4)
  end if
  
  !separate boundaries processors in different groups and communicators
  allocate(incletaN(1:dims(1)+1))
  !indices of processors in the topology that have real boundaries
  !include original master, to broadcast properly
  incletaN(1)=0
  do i=0,dims(1)-1
    incletaN(i+1+1)=(i+1)*dims(2)-1
  end do
  
  
  !get the original group
  call mpi_comm_group(comm2d,allgroup,ierror)
  
  !create groups for each boundary
  if (incletaN(2)==0) then
    call mpi_group_incl(allgroup,dims(1),incletaN(2:dims(1)+1),groupetaN,ierror)
  else
    call mpi_group_incl(allgroup,dims(1)+1,incletaN(1:dims(1)+1),groupetaN,ierror)
  end if
  
!   print*,incletaN
  
  !create communicators for each group  
  call mpi_comm_create(comm2d,groupetaN,commetaN,ierror)

  
  !catch boundary conditions
  if (any(myrank==incletaN)) then
!     print*,myrank
    select case (CB_real(4))
      case(0)
	!send parameters
	call mpi_bcast(nt_etaNg1,1,mpi_integer,master,commetaN,ierror)
	call mpi_bcast(dt_etaNg1,1,mpi_double_precision,master,commetaN,ierror)
	call mpi_bcast(optetaNg1,1,mpi_integer,master,commetaN,ierror)
	call mpi_bcast(nt_etaNg2,1,mpi_integer,master,commetaN,ierror)
	call mpi_bcast(dt_etaNg2,1,mpi_double_precision,master,commetaN,ierror)
	call mpi_bcast(optetaNg2,1,mpi_integer,master,commetaN,ierror) 
	call mpi_bcast(flagetaN,1,mpi_integer,master,commetaN,ierror) 
	
	!allocate buffer
	if (myrank/= master) then
	  allocate(bufqetaNg1(4,nxi(1),nt_etaNg1+1),bufqetaNg2(4,nxi(1),nt_etaNg2+1))
	end if
	!allocate a piece of bc
	allocate(qetaNg1(4,Nbx,nt_etaNg1),qetaNg2(4,Nbx,nt_etaNg2))
	call mpi_comm_rank(commetaN,myrank_temp,ierror)
	
	!bcast buffer
	call mpi_bcast(bufqetaNg1,4*nxi(1)*(nt_etaNg1+1),mpi_double_precision,master,commetaN,ierror)
	call mpi_bcast(bufqetaNg2,4*nxi(1)*(nt_etaNg2+1),mpi_double_precision,master,commetaN,ierror)  
	
	!initialize values
	qetaNg1=bufqetaNg1(:,si:ei,:)
	qetaNg2=bufqetaNg2(:,si:ei,:)
      
	!delete buffer
	deallocate(bufqetaNg1,bufqetaNg2)
	
	!catch exception (unnecesary master-flag)
	if (CB(4)/=0) then	
	  flagetaN=0
	end if
    end select 
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
end subroutine bcast_bcs

subroutine bcast_coupling_bc_1(ntg1,dtg1,optg1,ntg2,dtg2,optg2,flag,comm)
  !two soubroutines: need parameters to allocate, cant allocate inside a subroutine
  use mpi
  implicit none
  integer::ntg1,optg1,ntg2,optg2
  integer :: flag
  real (kind=8) ::dtg1,dtg2
  integer::comm,master,ierror  
  master=0  	
  call mpi_bcast(ntg1,1,mpi_integer,master,comm,ierror)
  call mpi_bcast(dtg1,1,mpi_double_precision,master,comm,ierror)
  call mpi_bcast(optg1,1,mpi_integer,master,comm,ierror)
  call mpi_bcast(ntg2,1,mpi_integer,master,comm,ierror)
  call mpi_bcast(dtg2,1,mpi_double_precision,master,comm,ierror)
  call mpi_bcast(optg2,1,mpi_integer,master,comm,ierror) 
  call mpi_bcast(flag,1,mpi_integer,master,comm,ierror) 

end subroutine bcast_coupling_bc_1

subroutine bcast_coupling_bc_2(ntg1,ntg2,buffqg1,buffqg2,qg1,qg2,s,e,N,M,comm,myrank1) 
  use mpi
  implicit none
  integer::ntg1,ntg2
  real (kind=8), dimension(4,N,ntg1+1) ::qg1
  real (kind=8), dimension(4,N,ntg2+1) ::qg2
  real (kind=8), dimension(4,M,ntg1+1) ::buffqg1
  real (kind=8), dimension(4,M,ntg2+1) ::buffqg2
  integer myrank1,myrank2
  integer::s,e,N,M,comm,master,ierror
  master=0  
  call mpi_comm_rank(comm,myrank2,ierror)


  call mpi_bcast(buffqg1,4*M*(ntg1+1),mpi_double_precision,master,comm,ierror)
  call mpi_bcast(buffqg2,4*M*(ntg2+1),mpi_double_precision,master,comm,ierror)  
  qg1=buffqg1(:,s:e,:)
  qg2=buffqg2(:,s:e,:)
  
  call mpi_barrier(comm,ierror)
end subroutine bcast_coupling_bc_2