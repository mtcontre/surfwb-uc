subroutine decomp2d
  USE global_variables
  USE geometries
  USE MPI
  USE MPI_SURF
  USE time0
  
  implicit none
  
  real(kind=8), dimension(:,:), allocatable :: buf_x,buf_y,buf_z, buf_fric
  real(kind=8), dimension(:,:,:), allocatable :: buf_qold
  character(len=255) :: fname,ofmt,command,filename
  logical::dir_e
  integer :: i
  
  !mpi knows how to split comm_world better into a cart topo:  
  !obtain dims vector given ndim and nproc (let mpi decide it)
  CALL MPI_DIMS_CREATE(nproc,ndim,dims,ierror)
  !create the communicator comm2d
  CALL MPI_CART_CREATE(MPI_COMM_WORLD,ndim,dims,isperiodic,&
    reorder,comm2d,ierror)
  !obtain myrank2d
  CALL MPI_COMM_RANK(COMM2D,myrank2d,ierror)
  !get the topology_coords in the cart topo for this process
  CALL MPI_CART_GET(COMM2D,ndim,dims,isperiodic,topology_coords,ierror)  
  
  !get my neighbors
  !left/right & back/front as when looking to the matrix
  !increasing indices starting from the bottom left corner
  !ex:  
  !	  COORDS
  !	1,0	1,1
  !^xi	0,0	0,1  
  !>eta
  !corresponds to
  !	myrank2d
  !	2	3
  !	0	1
  
  CALL MPI_CART_SHIFT(COMM2D,0,shift,myback,myfront,ierror) 
  CALL MPI_CART_SHIFT(COMM2D,1,shift,myleft,myright,ierror) 
  
  !obtain starting and ending grid indices si,ei and sj,ej
  CALL DECOMP1D(Nbx,dims(1),topology_coords(1),si,ei)
  CALL DECOMP1D(Nby,dims(2),topology_coords(2),sj,ej)    
  
  !re distribute pieces of 2darrays
  old_nbx = Nbx
  old_nby = Nby
  Nbx = ei-si+1
  Nby = ej-sj+1
  
  !assign boundary conditions  
  if (myback /= MPI_PROC_NULL) then
    CB(1) = -1
  end if
  
  if (myfront /= MPI_PROC_NULL) then
    CB(2) = -1
  end if
  
  if (myleft /= MPI_PROC_NULL) then
    CB(3) = -1
  end if
  
  if (myright /= MPI_PROC_NULL) then
    CB(4) = -1
  end if
  
  
  !save in buffers
  allocate(buf_x(Nbx,Nby), buf_y(Nbx,Nby), buf_z(Nbx,Nby), &
    buf_qold(3,Nbx,Nby), buf_fric(Nbx,Nby))
  buf_x = x_global(si:ei,sj:ej)
  buf_y = y_global(si:ei,sj:ej)
  buf_z = z_global(si:ei,sj:ej)
  buf_qold = qold_global(:,si:ei, sj:ej)
  buf_fric = Mcoef
  
  !up to this point, only matrices in mpisurf_input_geom.f90, 
  !mpisurf_ input_ic.f90 and Mcoef defined in init.f90 are defined.
  !We need to update them to have only the values of the processor's grid
  
  !Forget old and define new local 2darrays
  deallocate(x_global, y_global, z_global, qold_global, qnew_global, q0_global, &
     qreal_global, V_global, C_global, VC,S1_global,S2_global, Mcoef)
  allocate(x_global(Nbx,Nby), y_global(Nbx,Nby), z_global(Nbx,Nby), &
    qold_global(3,Nbx,Nby), qnew_global(3,Nbx,Nby), q0_global(3,Nbx,Nby), &
     qreal_global(3,Nbx,Nby), V_global(Nbx,Nby), C_global(Nbx,Nby), &
     VC(Nbx,Nby),S1_global(Nbx,Nby),S2_global(Nbx,Nby), Mcoef(Nbx,Nby))
  !im ommiting GA variables defined in input_ic.f90 ...(?)
  
  x_global = buf_x
  y_global = buf_y
  z_global = buf_z
  qold_global = buf_qold     
  Mcoef = buf_fric
  
  !forget buffers
  deallocate(buf_x, buf_y, buf_z, buf_qold, Mcoef)
  
  !now save topology data for posterior output reconstruction
  if (myrank==0) then    
    !check if outdir/grids exists
    fname=trim(outdir)//'/grids/'
    inquire(file=fname,exist=dir_e)    
    if( .not. dir_e) then
      command='mkdir '//trim(fname)
      call system(trim(command))
    end if   

    !write gridproperties.dat
    fname=trim(outdir)//'/grids/gridproperties.dat'
    open(unit=50,file=fname)
    write(unit=50,fmt=*)'dit', dit
    write(unit=50,fmt='("nproc ",I3.3)') nproc  
    write(unit=50,fmt='("dims ",I3.3,X,I3.3)') dims(1),dims(2)
    write(unit=50,fmt='("nxi ", I3.3)'),old_nbx
    write(unit=50,fmt='("neta ", I3.3)')old_nby 
    write(unit=50,fmt=*) trim(batiname(1))
    write(unit=50,fmt=*) trim(batiname(2))
    write(unit=50,fmt=*) trim(batiname(3))
    
    do i=1,3
      command='cp '//trim(batiname(i))//' '//trim(outdir)//'/.'
      call system(command)
    end do   
  end if
  
  !write this grid
  write(filename,'("/grids/grid",I3.3,"_",I3.3,".dat")') topology_coords(1),topology_coords(2)
  filename=trim(outdir)//trim(adjustl(filename))
  open(unit=myrank+100,file=filename)
  write(unit=myrank+100,fmt='( I3.3, "  Nbx" )') Nbx
  write(unit=myrank+100,fmt='( I3.3, "  Nby" )') Nby  
  write(unit=myrank+100,fmt='( I3.3, "  si" )') si
  write(unit=myrank+100,fmt='( I3.3, "  ei" )') ei
  write(unit=myrank+100,fmt='( I3.3, "  sj" )') sj
  write(unit=myrank+100,fmt='( I3.3, "  ej" )') ej
  write(unit=myrank+100,fmt='( I3.3, "  coord(1)" )') topology_coords(1)
  write(unit=myrank+100,fmt='( I3.3, "  coord(2)" )') topology_coords(2)
  write(unit=myrank+100,fmt='( I3.3, "  dims(1)" )') dims(1)
  write(unit=myrank+100,fmt='( I3.3, "  dims(2)" )') dims(2)  
  write(unit=myrank+100,fmt='( I3.3, "  rank2d" )') myrank2d
  write(unit=myrank+100,fmt='( I4.3, "  left" )') myleft
  write(unit=myrank+100,fmt='( I4.3, "  right" )') myright
  write(unit=myrank+100,fmt='( I4.3, "  back" )') myback
  write(unit=myrank+100,fmt='( I4.3, "  front" )') myfront
  close(unit=myrank+100)
end subroutine

subroutine decomp1d(n,p,rank,s,e)
  !recibe el numero total n y el numero de bloques p
  !rank es el indice del que queremos saber
  !los numeros s y e
  !s es el indice inferior
  !e es el indice superior
  !la forma mas facil es ver que pasa con s_rank, y e=s_{rank+1}-1
  !para los diferentes casos: rank<deficit o rank>=deficit
  !si rank<deficit s=s+r, sino s=s+deficit....esa es la diferencia
  !la idea es pasar el deficit a los primeros procesadores y asi balancear la cosa
  implicit none
  integer::n,p,rank,s,e
  integer::nlocal,deficit
  nlocal=n/p
  s=rank*nlocal+1
  deficit=mod(n,p)
  if (rank<deficit) then
    s=s+rank
    nlocal=nlocal+1
  else
    s=s+deficit
  end if
  e=s+nlocal-1
  if (e>n .or. rank==p-1) then
    e=n
  end if     
end subroutine