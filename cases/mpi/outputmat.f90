subroutine outputmat
  !prints everything and q0 to files
  use global_variables
  use geometries
  use timeseries
  
  implicit none
  integer::i,j
  character(len=255) ::filename,filenameT,ofmt,fout_param
  logical::dir_exists,t_exists,param_exist
  logical::to_out=.false.
  character(len=8)::charit
  if (dit==-1) then
    to_out = print_out
    print_out=.False.
  else
    to_out=mod(it,dit)==0  
    nitout = it
  end if
  if (to_out.or.treal==tinit) then
    filenameT=trim(outdir)//'/time.dat'
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
  
  !allocate at first time step
  if (treal==tinit) then
    qnew_global(:,:,:)=qold_global(:,:,:)
  end if
  
  !dimensionalize to save results
  qreal_global(1,:,:)=qnew_global(1,:,:)*H
  qreal_global(2,:,:)=qnew_global(2,:,:)*U
  qreal_global(3,:,:)=qnew_global(3,:,:)*U  
  
  !write results to file
  if (to_out.or.treal==tinit) then    
    write(charit,'(I8.8)') nitout
      filename=trim(outdir)//'/SOL2D.'//trim(charit)//'.dat'
!     filename='results/SOL2D.'//trim(charit)//'.dat'
!     filename=trim(filename)//'.dat'
    open(unit=100,file=filename)
    do j=1,Nby
      do i=1,Nbx
	write(unit=100,fmt='(6(E20.10,2X))') &
	    x_global(i,j),y_global(i,j),z_global(i,j),qold_global(1,i,j), &
	    qold_global(2,i,j),qold_global(3,i,j)
      end do
    end do
    close(unit=100)
    
    !write parameters
    fout_param=trim(outdir)//'/param.dat'
    inquire(file=fout_param, exist=param_exist)
    if (param_exist) then
      open (25,file=fout_param, status='replace',action='write')
    else
      open(25,file=fout_param, status='new', action='write')
    end if
    
    write(25,'(6(I5,2x))') caso, Nbx, Nby, Nts,nitout, dit
    close(25)
    
    call system('gzip -f '//filename)
  end if
  

end subroutine