subroutine outputmat_par
  !prints everything and q0 to files
  use global_variables
  use geometries
  use mpi
  use mpi_surf
  use multigrid_surf
  implicit none
  integer::i,j
  character(len=255) ::filename,filenameT,ofmt
  logical::dir_exists,t_exists
  logical::to_out=.false.

  !determine if should save output, and restart flag
  if (dit==-1) then
    to_out = print_out
    print_out=.False.
  else
    to_out=mod(it,dit)==0  
    nitout = it
  end if  
  
  !save time steps list
  if (myrank==master) then
!     filenameT='results/time.dat'
    if (to_out.or.treal==tinit) then
      write(filenameT,'("/timeP",I3.3,".dat")')nproc
      filenameT=trim(outdir)//trim(filenameT)
      inquire(file=filenameT,exist=t_exists)
      if (.not. t_exists) then
	open(1,file=filenameT,status='new',action='write')
      else
	if (treal==tinit) then
	  open(1,file=filenameT,status='replace',action='write')
	else
	  open(1,file=filenameT,status='old',action='write',position='append')
	end if
      end if
      write(unit=1,fmt='(E20.10)') treal
      close(1)  
    end if
  end if
  
  !allocate variables at first time step
  
  if (treal==tinit) then
    allocate(qnew_global(3,Nbx,Nby),qreal_global(3,Nbx,Nby))
    qnew_global(:,:,:)=qold_global(:,:,:)
  end if
  
  !dimensionalize to save results
  qreal_global(1,:,:)=qnew_global(1,:,:)*H
  qreal_global(2,:,:)=qnew_global(2,:,:)*U
  qreal_global(3,:,:)=qnew_global(3,:,:)*U  
  
  !write results to file
  if (to_out.or.treal==tinit)then    
    write(filename,'("/P",I3.3,"frame",I8.8,".",I3.3,"_",I3.3,".dat")') nproc,nitout,coords(1),coords(2)
    filename=trim(adjustl(outdir))//trim(filename)
    open(unit=myrank+100,file=filename)
    do j=1,Nby
      do i=1,Nbx
	write(unit=myrank+100,fmt='(6(E20.10,2X))') qold_global(1,i,j), &
	  qold_global(1,i,j)+z_global(i,j),qold_global(2,i,j),qold_global(3,i,j)
      end do
    end do
    close(unit=myrank+100)
  end if
end subroutine