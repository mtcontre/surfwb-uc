subroutine exchange_2d(qt)
  use mpi
  use mpi_surf
  use global_variables
  integer,dimension(mpi_status_size)::status
  real(kind=8),dimension(3,Nbx+4,Nby+4) :: qt
!   print*,myleft,myrank2d,myright
  !the real domain for each grid is in range (3:Nbx+2,3:Nby+2) 
  !the other rows represent boundary ghost cells

  !send xi_(n+1) to my front, recv xi_1 from my back .....the tag coordinates them
  call mpi_sendrecv(qt(:,Nbx+1,3:Nby+2),3*Nby,mpi_double_precision,myfront,100, &	
	qt(:,1,3:Nby+2),3*Nby,mpi_double_precision,myback,100,comm2d,status,ierror)
  
  !send xi_(n+2) to my front, recv xi_2 from myback
  call mpi_sendrecv(qt(:,Nbx+2,3:Nby+2),3*Nby,mpi_double_precision,myfront,101, &	
	qt(:,2,3:Nby+2),3*Nby,mpi_double_precision,myback,101,comm2d,status,ierror)
  
  !send xi_3 to myback,  recv xi_(n+3) from my front
  call mpi_sendrecv(qt(:,3,3:Nby+2),3*Nby,mpi_double_precision,myback,102, &	
	qt(:,Nbx+3,3:Nby+2),3*Nby,mpi_double_precision,myfront,102,comm2d,status,ierror)
  
  !send xi_4 to myback, recv xi_(n+4) from my front
  call mpi_sendrecv(qt(:,4,3:Nby+2),3*Nby,mpi_double_precision,myback,103, &	
	qt(:,Nbx+4,3:Nby+2),3*Nby,mpi_double_precision,myfront,103,comm2d,status,ierror)
	
	
  !send eta_(n+1) to my right, recv eta_1 from myleft
  call mpi_sendrecv(qt(:,3:Nbx+2,Nby+1),3*Nbx,mpi_double_precision,myright,200, &	
	qt(:,3:Nbx+2,1),3*Nbx,mpi_double_precision,myleft,200,comm2d,status,ierror)
	
  !send eta_(n+2) to my right, recv eta_2 from myleft
  call mpi_sendrecv(qt(:,3:Nbx+2,Nby+2),3*Nbx,mpi_double_precision,myright,201, &	
	qt(:,3:Nbx+2,2),3*Nbx,mpi_double_precision,myleft,201,comm2d,status,ierror)
	
  !send eta_3 to myleft, recv eta_(n+3) from my right
  call mpi_sendrecv(qt(:,3:Nbx+2,3),3*Nbx,mpi_double_precision,myleft,202, &	
	qt(:,3:Nbx+2,Nby+3),3*Nbx,mpi_double_precision,myright,202,comm2d,status,ierror)
  !send eta_4 to myleft, recv eta_(n+4)from my right
  call mpi_sendrecv(qt(:,3:Nbx+2,4),3*Nbx,mpi_double_precision,myleft,203, &	
	qt(:,3:Nbx+2,Nby+4),3*Nbx,mpi_double_precision,myright,203,comm2d,status,ierror)
end subroutine
! 
! call mpi_sendrecv(qt(:,Nbx+1,3:Nby+2),3*Nby,mpi_double_precision,myright,100, &	
! 	qt(:,1,3:Nby+2),3*Nby,mpi_double_precision,myleft,100,comm2d,status,ierror)
!   
!   call mpi_sendrecv(qt(:,Nbx+2,3:Nby+2),3*Nby,mpi_double_precision,myright,101, &	
! 	qt(:,2,3:Nby+2),3*Nby,mpi_double_precision,myleft,101,comm2d,status,ierror)
!   
!   !send xi 1, xi 2, to myleft. receive xiN,xiN+1 from my right
!   call mpi_sendrecv(qt(:,1,3:Nby+2),3*Nby,mpi_double_precision,myleft,102, &	
! 	qt(:,Nbx+1,3:Nby+2),3*Nby,mpi_double_precision,myright,102,comm2d,status,ierror)
!   
!   call mpi_sendrecv(qt(:,2,3:Nby+2),3*Nby,mpi_double_precision,myleft,103, &	
! 	qt(:,Nbx+2,3:Nby+2),3*Nby,mpi_double_precision,myright,103,comm2d,status,ierror)
! 	
! 	
!   !send eta N-1, eta N to myfront, receive eta 1,eta  2 from myback
!   call mpi_sendrecv(qt(:,3:Nbx+2,Nby+1),3*Nbx,mpi_double_precision,myfront,200, &	
! 	qt(:,3:Nbx+2,1),3*Nbx,mpi_double_precision,myback,200,comm2d,status,ierror)
! 	
!   call mpi_sendrecv(qt(:,3:Nbx+2,Nby+2),3*Nbx,mpi_double_precision,myfront,201, &	
! 	qt(:,3:Nbx+2,2),3*Nbx,mpi_double_precision,myback,201,comm2d,status,ierror)
! 	
!   !send eta 1, eta 2 to myback, receive eta N-1,N from my front
!   call mpi_sendrecv(qt(:,3:Nbx+2,1),3*Nbx,mpi_double_precision,myback,202, &	
! 	qt(:,3:Nbx+2,Nby+1),3*Nbx,mpi_double_precision,myfront,202,comm2d,status,ierror)
! 	
!   call mpi_sendrecv(qt(:,3:Nbx+2,2),3*Nbx,mpi_double_precision,myback,203, &	
! 	qt(:,3:Nbx+2,Nby+2),3*Nbx,mpi_double_precision,myfront,203,comm2d,status,ierror)