subroutine outputmat_par
  !prints everything and q0 to files
  use global_variables
  use geometries
  use mpi
  use mpi_surf
  implicit none
  integer::i,j
  character(len=255) ::filename,filenameT,ofmt,command
  logical::dir_exists,t_exists, to_out=.false.
  
  !create outdir if neccesary
  inquire(file=trim(outdir),exist=dir_exists)
  if (.not. dir_exists) then
    command = 'mkdir '//trim(outdir)
    call system(trim(command))
  end if
  
  !check if need to print
  if (dit==-1) then
    to_out = print_out
    print_out=.False.
  else
    to_out=mod(it,dit)==0  
    nitout = it
  end if
  
  !master fills time file
  if (myrank==master) then
    if (mod(it,dit)==0) then
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
  
  !dimensionalize to save results
  qreal_global(1,:,:)=qnew_global(1,:,:)*H
  qreal_global(2,:,:)=qnew_global(2,:,:)*U
  qreal_global(3,:,:)=qnew_global(3,:,:)*U  
  
  if (to_out) then    
    write(filename,'("/P",I3.3,"frame",I8.8,".",I3.3,"_",I3.3,".dat")') nproc,it,coords(1),coords(2)
    filename=trim(adjustl(outdir))//trim(filename)
    open(unit=myrank+100,file=filename)
    
    do i=1,Nbx
      do j=1,Nby
	write(unit=myrank+100,fmt='(3(E20.10,2X))') qold_global(1,i,j), qold_global(2,i,j),qold_global(3,i,j)
      end do
    end do
    close(unit=myrank+100)
  end if
end subroutine