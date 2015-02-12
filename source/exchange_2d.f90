subroutine exchange_2d(qt,xit,etat,zt)
  use mpi
  use mpi_surf
  use global_variables,only:Nbx,Nby
  implicit none
  real(kind=8), dimension(3,Nbx+4,Nby+4), intent(inout) :: qt
  real(kind=8), dimension(2,Nbx+4,Nby+4), intent(inout) ::xit,etat
  real(kind=8), dimension(Nbx+4,Nby+4), intent(inout) ::zt

  integer, dimension(mpi_status_size)::status
  real(kind=8), dimension(2*(Nby+4)) :: buff,bufb,fbuf, bbuf
  real(kind=8), dimension(2*(Nbx+4)) :: bufl, bufr, lbuf, rbuf
  integer :: nvar
  
  !-------exchange q------
  
  do nvar=1,3
    !pack
    bufb(1:(Nby+4))	= qt(nvar,3,:)			!back
    bufb((Nby+5):2*(Nby+4)) = qt(nvar,4,:)
    buff(1:(Nby+4))	= qt(nvar,Nbx+1,:)		!front
    buff((Nby+5):2*(Nby+4)) = qt(nvar,Nbx+2, :)
    
    bufr(1:(Nbx+4))	= qt(nvar,:,Nby+1)		!right
    bufr((Nbx+5):2*(Nbx+4)) = qt(nvar,:,Nby+2)
    bufl(1:(Nbx+4))	= qt(nvar,:,3)			!left
    bufl((Nbx+5):2*(Nbx+4)) = qt(nvar,:,4)
    
    !exchange packed data
    call mpi_sendrecv(buff, 2*(Nby+4), mpi_double_precision, myfront, 1+10*(nvar-1), &
		      bbuf, 2*(Nby+4), mpi_double_precision, myback, 1+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufb, 2*(Nby+4), mpi_double_precision, myback, 2+10*(nvar-1), &
		      fbuf, 2*(Nby+4), mpi_double_precision, myfront, 2+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufr, 2*(Nbx+4), mpi_double_precision, myright, 3+10*(nvar-1), &
		      lbuf, 2*(Nbx+4), mpi_double_precision, myleft, 3+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufl, 2*(Nbx+4), mpi_double_precision, myleft, 4+10*(nvar-1), &
		      rbuf, 2*(Nbx+4), mpi_double_precision, myright, 4+10*(nvar-1), &
		      comm2d, status, ierror)
    !unpack
    if (myfront /= mpi_proc_null) then
      qt(nvar,Nbx+3,:) = fbuf(1:(Nby+4))
      qt(nvar,Nbx+4,:) = fbuf((Nby+5):2*(Nby+4))
    end if
    
    if (myback /= mpi_proc_null) then
      qt(nvar,1,:) = bbuf(1:(Nby+4))
      qt(nvar,2,:) = bbuf((Nby+5):2*(Nby+4))
    end if
    
    if (myright /= mpi_proc_null) then
      qt(nvar,:,Nby+3) = rbuf(1:(Nbx+4))
      qt(nvar,:,Nby+4) = rbuf((Nbx+5):2*(Nbx+4))
    end if
    
    if (myleft /= mpi_proc_null) then
      qt(nvar,:,1) = lbuf(1:(Nbx+4))
      qt(nvar,:,2) = lbuf((Nbx+5):2*(Nbx+4))
    end if  
  end do

  !-------exchange xit------
  do nvar=1,2
    !pack
    bufb(1:(Nby+4))	= xit(nvar,3,:)			!back
    bufb((Nby+5):2*(Nby+4)) = xit(nvar,4,:)
    buff(1:(Nby+4))	= xit(nvar,Nbx+1,:)		!front
    buff((Nby+5):2*(Nby+4)) = xit(nvar,Nbx+2, :)
    
    bufr(1:(Nbx+4))	= xit(nvar,:,Nby+1)		!right
    bufr((Nbx+5):2*(Nbx+4)) = xit(nvar,:,Nby+2)
    bufl(1:(Nbx+4))	= xit(nvar,:,3)		!left
    bufl((Nbx+5):2*(Nbx+4)) = xit(nvar,:,4)
    
    !exchange packed data
    call mpi_sendrecv(buff, 2*(Nby+4), mpi_double_precision, myfront, 31+10*(nvar-1), &
		      bbuf, 2*(Nby+4), mpi_double_precision, myback, 31+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufb, 2*(Nby+4), mpi_double_precision, myback, 32+10*(nvar-1), &
		      fbuf, 2*(Nby+4), mpi_double_precision, myfront, 32+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufr, 2*(Nbx+4), mpi_double_precision, myright, 33+10*(nvar-1), &
		      lbuf, 2*(Nbx+4), mpi_double_precision, myleft, 33+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufl, 2*(Nbx+4), mpi_double_precision, myleft, 34+10*(nvar-1), &
		      rbuf, 2*(Nbx+4), mpi_double_precision, myright, 34+10*(nvar-1), &
		      comm2d, status, ierror)
    !unpack
    if (myfront /= mpi_proc_null) then
      xit(nvar,Nbx+3,:) = fbuf(1:(Nby+4))
      xit(nvar,Nbx+4,:) = fbuf((Nby+5):2*(Nby+4))
    end if
    
    if (myback /= mpi_proc_null) then
      xit(nvar,1,:) = bbuf(1:(Nby+4))
      xit(nvar,2,:) = bbuf((Nby+5):2*(Nby+4))
    end if
    
    if (myright /= mpi_proc_null) then
      xit(nvar,:,Nby+3) = rbuf(1:(Nbx+4))
      xit(nvar,:,Nby+4) = rbuf((Nbx+5):2*(Nbx+4))
    end if
    
    if (myleft /= mpi_proc_null) then
      xit(nvar,:,1) = lbuf(1:(Nbx+4))
      xit(nvar,:,2) = lbuf((Nbx+5):2*(Nbx+4))
    end if  
  end do
  
  
  
  
  !-------exchange etat------
  do nvar=1,2
    !pack
    bufb(1:(Nby+4))	= etat(nvar,3,:)		!back
    bufb((Nby+5):2*(Nby+4)) = etat(nvar,4,:)
    buff(1:(Nby+4))	= etat(nvar,Nbx+1,:)		!front
    buff((Nby+5):2*(Nby+4)) = etat(nvar,Nbx+2, :)
    
    bufr(1:(Nbx+4))	= etat(nvar,:,Nby+1)		!right
    bufr((Nbx+5):2*(Nbx+4)) = etat(nvar,:,Nby+2)
    bufl(1:(Nbx+4))	= etat(nvar,:,3)		!left
    bufl((Nbx+5):2*(Nbx+4)) = etat(nvar,:,4)
    
    !exchange packed data
    call mpi_sendrecv(buff, 2*(Nby+4), mpi_double_precision, myfront, 51+10*(nvar-1), &
		      bbuf, 2*(Nby+4), mpi_double_precision, myback, 51+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufb, 2*(Nby+4), mpi_double_precision, myback, 52+10*(nvar-1), &
		      fbuf, 2*(Nby+4), mpi_double_precision, myfront, 52+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufr, 2*(Nbx+4), mpi_double_precision, myright, 53+10*(nvar-1), &
		      lbuf, 2*(Nbx+4), mpi_double_precision, myleft, 53+10*(nvar-1), &
		      comm2d, status, ierror)
    call mpi_sendrecv(bufl, 2*(Nbx+4), mpi_double_precision, myleft, 54+10*(nvar-1), &
		      rbuf, 2*(Nbx+4), mpi_double_precision, myright, 54+10*(nvar-1), &
		      comm2d, status, ierror)
    !unpack
    if (myfront /= mpi_proc_null) then
      etat(nvar,Nbx+3,:) = fbuf(1:(Nby+4))
      etat(nvar,Nbx+4,:) = fbuf((Nby+5):2*(Nby+4))
    end if
    
    if (myback /= mpi_proc_null) then
      etat(nvar,1,:) = bbuf(1:(Nby+4))
      etat(nvar,2,:) = bbuf((Nby+5):2*(Nby+4))
    end if
    
    if (myright /= mpi_proc_null) then
      etat(nvar,:,Nby+3) = rbuf(1:(Nbx+4))
      etat(nvar,:,Nby+4) = rbuf((Nbx+5):2*(Nbx+4))
    end if
    
    if (myleft /= mpi_proc_null) then
      etat(nvar,:,1) = lbuf(1:(Nbx+4))
      etat(nvar,:,2) = lbuf((Nbx+5):2*(Nbx+4))
    end if  
  end do
  
  
  
  !-------exchange zt------
  !--------zt could be constant through the simmulation
  !-------there seems to be no need to communicate so many times
  !-------unless having variable zt (?)
  
  !pack
  bufb(1:(Nby+4))	= zt(3,:)		!back
  bufb((Nby+5):2*(Nby+4)) = zt(4,:)
  buff(1:(Nby+4))	= zt(Nbx+1,:)		!front
  buff((Nby+5):2*(Nby+4)) = zt(Nbx+2, :)
  
  bufr(1:(Nbx+4))	= zt(:,Nby+1)		!right
  bufr((Nbx+5):2*(Nbx+4)) = zt(:,Nby+2)
  bufl(1:(Nbx+4))	= zt(:,3)		!left
  bufl((Nbx+5):2*(Nbx+4)) = zt(:,4)
  
  !exchange packed data
  call mpi_sendrecv(buff, 2*(Nby+4), mpi_double_precision, myfront, 71+10*(nvar-1), &
		    bbuf, 2*(Nby+4), mpi_double_precision, myback, 71+10*(nvar-1), &
		    comm2d, status, ierror)
  call mpi_sendrecv(bufb, 2*(Nby+4), mpi_double_precision, myback, 72+10*(nvar-1), &
		    fbuf, 2*(Nby+4), mpi_double_precision, myfront, 72+10*(nvar-1), &
		    comm2d, status, ierror)
  call mpi_sendrecv(bufr, 2*(Nbx+4), mpi_double_precision, myright, 73+10*(nvar-1), &
		    lbuf, 2*(Nbx+4), mpi_double_precision, myleft, 73+10*(nvar-1), &
		    comm2d, status, ierror)
  call mpi_sendrecv(bufl, 2*(Nbx+4), mpi_double_precision, myleft, 74+10*(nvar-1), &
		    rbuf, 2*(Nbx+4), mpi_double_precision, myright, 74+10*(nvar-1), &
		    comm2d, status, ierror)
  !unpack
  if (myfront /= mpi_proc_null) then
    zt(Nbx+3,:) = fbuf(1:(Nby+4))
    zt(Nbx+4,:) = fbuf((Nby+5):2*(Nby+4))
  end if
  
  if (myback /= mpi_proc_null) then
    zt(1,:) = bbuf(1:(Nby+4))
    zt(2,:) = bbuf((Nby+5):2*(Nby+4))
  end if
  
  if (myright /= mpi_proc_null) then
    zt(:,Nby+3) = rbuf(1:(Nbx+4))
    zt(:,Nby+4) = rbuf((Nbx+5):2*(Nbx+4))
  end if
  
  if (myleft /= mpi_proc_null) then
    zt(:,1) = lbuf(1:(Nbx+4))
    zt(:,2) = lbuf((Nbx+5):2*(Nbx+4))
  end if  
  
  
!   print*,'Communication done on proc',myrank
  
end subroutine
