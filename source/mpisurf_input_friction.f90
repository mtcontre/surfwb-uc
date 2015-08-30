subroutine input_friction
  use mpi_surf
  use multigrid_surf
  use global_variables
  implicit none  
  integer::i,j,level
  
  allocate(buffMcoef(ngrids))
  do level=1,ngrids
    allocate(buffMcoef(level)%bMcoef(nxi(level),neta(level)))
    if (fopt==0) then
      do i=1,nxi(level)
	do j=1,neta(level)
	  buffMcoef(level)%bMcoef(i,j)=0.0D0
	end do
      end do
    else
      if (fM==1) then!single value
	do i=1,nxi(level)
	  do j=1,neta(level)
	    buffMcoef(level)%bMcoef(i,j)=Coef
	  end do
	end do
      else if (fM==2) then!matrix
	open	(unit=99, file ='data/friction_run31.dat')
	read(99,*) ((buffMCoef(level)%bMcoef(i,j),i=1,nxi(level)),j=1,neta(level))
	close(unit=99)
      end if
    end if
  end do  
end subroutine input_friction
