subroutine outputmat_par
  !prints everything and q0 to files
  use global_variables
  use geometries
  use mpi
  use mpi_surf
  use multigrid_surf
  implicit none
  integer::i,j
  character(len=255) ::filename,ofmt
  logical::dir_e
  if (myrank==master) then
    open(unit=50,file='results/gridproperties.dat')
    write(unit=50,fmt='("nproc=",I3.3)') nproc  
    write(unit=50,fmt='("dims=",I3.3,X,I3.3)') dims(1),dims(2)
    write(unit=50,fmt='("ngrids=",I3.3)') ngrids
    write(ofmt,'("(",I3,"(A5,I4.4,X))")') ngrids
    write(unit=50,fmt=ofmt)'nxi =',nxi
    write(unit=50,fmt=ofmt)'neta=',neta    
    close(unit=50)    
    inquire(file='results/grids/.',exist=dir_e)    
    if( .not. dir_e) then
      call system('mkdir results/grids')
    end if    
  end if
  
  !allocate at first time step
  if (treal==tinit) then
    allocate(qnew_global(3,Nbx,Nby),qreal_global(3,Nbx,Nby))
    qnew_global(:,:,:)=qold_global(:,:,:)
  end if
  
  !dimensionalize to save results
  qreal_global(1,:,:)=qnew_global(1,:,:)*H
  qreal_global(2,:,:)=qnew_global(2,:,:)*U
  qreal_global(3,:,:)=qnew_global(3,:,:)*U  
  
  if (mod(it,dit)==0) then
    write(filename,'("results/frame",I8.8,".",I3.3,"_",I3.3,".dat")') it,coords(1),coords(2)
    open(unit=myrank+100,file=filename)
    do i=1,Nbx
      do j=1,Nbx
	write(unit=myrank+100,fmt='(6(E15.5,2X))') qold_global(1,i,j),qold_global(1,i,j),qold_global(1,i,j)
      end do
    end do
    close(unit=myrank+100)
  end if
end subroutine