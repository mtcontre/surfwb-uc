subroutine exchange_2d(qt,xit,etat,zt)
  use mpi
  use mpi_surf
  use global_variables,only:Nbx,Nby
  integer,dimension(mpi_status_size)::status
  real(kind=8),dimension(3,Nbx+4,Nby+4) :: qt
  real(kind=8),dimension(2,Nbx+4,Nby+4) ::xit,etat
  real(kind=8),dimension(Nbx+4,Nby+4) ::zt
  !the real domain for each grid is in range (3:Nbx+2,3:Nby+2) 
  !the other rows represent boundary ghost cells  
  
  !sendrecv q
  call exchng_var_qxieta(qt,3)
  
  !sendrecv xi
  call exchng_var_qxieta(xit,2)
  
! !   sendrecv eta
  call exchng_var_qxieta(etat,2)
! !   
! !   sendrecv zt
  call exchng_var_z(zt)
  
end subroutine

subroutine  exchng_var_qxieta(var,dim1)
  use mpi
  use mpi_surf
  use global_variables,only:Nbx,Nby
  integer::dim1
  integer,dimension(mpi_status_size)::status
  real(kind=8),dimension(dim1,Nbx+4,Nby+4) :: var
  !the real domain for each grid is in range (3:Nbx+2,3:Nby+2) 
  !the other rows represent boundary ghost cells

  
  !***********sendrecv  q*************
  !send xi_(n+1) to my front, recv xi_1 from my back .....the tag coordinates them
  call mpi_sendrecv(var(:,Nbx+1,3:Nby+2),dim1*Nby,mpi_double_precision,myfront,100, &	
	var(:,1,3:Nby+2),dim1*Nby,mpi_double_precision,myback,100,comm2d,status,ierror)
  
  !send xi_(n+2) to my front, recv xi_2 from myback
  call mpi_sendrecv(var(:,Nbx+2,3:Nby+2),dim1*Nby,mpi_double_precision,myfront,101, &	
	var(:,2,3:Nby+2),dim1*Nby,mpi_double_precision,myback,101,comm2d,status,ierror)
  
  !send xi_3 to myback,  recv xi_(n+3) from my front
  call mpi_sendrecv(var(:,3,3:Nby+2),dim1*Nby,mpi_double_precision,myback,102, &	
	var(:,Nbx+3,3:Nby+2),dim1*Nby,mpi_double_precision,myfront,102,comm2d,status,ierror)
  
  !send xi_4 to myback, recv xi_(n+4) from my front
  call mpi_sendrecv(var(:,4,3:Nby+2),dim1*Nby,mpi_double_precision,myback,103, &	
	var(:,Nbx+4,3:Nby+2),dim1*Nby,mpi_double_precision,myfront,103,comm2d,status,ierror)
	
	
  !send eta_(n+1) to my right, recv eta_1 from myleft
  call mpi_sendrecv(var(:,3:Nbx+2,Nby+1),dim1*Nbx,mpi_double_precision,myright,200, &	
	var(:,3:Nbx+2,1),dim1*Nbx,mpi_double_precision,myleft,200,comm2d,status,ierror)
	
  !send eta_(n+2) to my right, recv eta_2 from myleft
  call mpi_sendrecv(var(:,3:Nbx+2,Nby+2),dim1*Nbx,mpi_double_precision,myright,201, &	
	var(:,3:Nbx+2,2),dim1*Nbx,mpi_double_precision,myleft,201,comm2d,status,ierror)
	
  !send eta_3 to myleft, recv eta_(n+3) from my right
  call mpi_sendrecv(var(:,3:Nbx+2,3),dim1*Nbx,mpi_double_precision,myleft,202, &	
	var(:,3:Nbx+2,Nby+3),dim1*Nbx,mpi_double_precision,myright,202,comm2d,status,ierror)
  !send eta_4 to myleft, recv eta_(n+4)from my right
  call mpi_sendrecv(var(:,3:Nbx+2,4),dim1*Nbx,mpi_double_precision,myleft,203, &	
	var(:,3:Nbx+2,Nby+4),dim1*Nbx,mpi_double_precision,myright,203,comm2d,status,ierror)
  
end subroutine

subroutine  exchng_var_z(var)
  use mpi
  use mpi_surf
  use global_variables,only:Nbx,Nby
!   integer::dim1
  integer,dimension(mpi_status_size)::status
  real(kind=8),dimension(Nbx+4,Nby+4) :: var
  !the real domain for each grid is in range (3:Nbx+2,3:Nby+2) 
  !the other rows represent boundary ghost cells
  

  !send xi_(n+1) to my front, recv xi_1 from my back .....the tag coordinates them
  call mpi_sendrecv(var(Nbx+1,3:Nby+2),Nby,mpi_double_precision,myfront,100, &	
	var(1,3:Nby+2),Nby,mpi_double_precision,myback,100,comm2d,status,ierror)
  
  !send xi_(n+2) to my front, recv xi_2 from myback
  call mpi_sendrecv(var(Nbx+2,3:Nby+2),Nby,mpi_double_precision,myfront,101, &	
	var(2,3:Nby+2),Nby,mpi_double_precision,myback,101,comm2d,status,ierror)
  
  !send xi_3 to myback,  recv xi_(n+3) from my front
  call mpi_sendrecv(var(3,3:Nby+2),Nby,mpi_double_precision,myback,102, &	
	var(Nbx+3,3:Nby+2),Nby,mpi_double_precision,myfront,102,comm2d,status,ierror)
  
  !send xi_4 to myback, recv xi_(n+4) from my front
  call mpi_sendrecv(var(4,3:Nby+2),Nby,mpi_double_precision,myback,103, &	
	var(Nbx+4,3:Nby+2),Nby,mpi_double_precision,myfront,103,comm2d,status,ierror)
	
	
  !send eta_(n+1) to my right, recv eta_1 from myleft
  call mpi_sendrecv(var(3:Nbx+2,Nby+1),Nbx,mpi_double_precision,myright,200, &	
	var(3:Nbx+2,1),Nbx,mpi_double_precision,myleft,200,comm2d,status,ierror)
	
  !send eta_(n+2) to my right, recv eta_2 from myleft
  call mpi_sendrecv(var(3:Nbx+2,Nby+2),Nbx,mpi_double_precision,myright,201, &	
	var(3:Nbx+2,2),Nbx,mpi_double_precision,myleft,201,comm2d,status,ierror)
	
  !send eta_3 to myleft, recv eta_(n+3) from my right
  call mpi_sendrecv(var(3:Nbx+2,3),Nbx,mpi_double_precision,myleft,202, &	
	var(3:Nbx+2,Nby+3),Nbx,mpi_double_precision,myright,202,comm2d,status,ierror)
  !send eta_4 to myleft, recv eta_(n+4)from my right
  call mpi_sendrecv(var(3:Nbx+2,4),Nbx,mpi_double_precision,myleft,203, &	
	var(3:Nbx+2,Nby+4),Nbx,mpi_double_precision,myright,203,comm2d,status,ierror)
  
end subroutine
