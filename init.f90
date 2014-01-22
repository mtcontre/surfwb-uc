!Rutina que inicia las variables,lee datos de entrada, lee batimetria, condiciones iniciales, condiciones de borde y asigna estos datos a variables
!Llama a tranformacion de coordenadas
SUBROUTINE init
  use mpi
  use mpi_surf
  use global_variables
  
  use custombc
  use global_variables
  use geometries
  use coords
  implicit none

  integer :: i,j 
  real (kind=8)::U1,U2

  !some relevant parameters
  g=9.812D0
  hmin=1.0E-7
  it=0.0D0
  dt=0.0D0
  t=0.0D0!0.0D0
  treal=0.0D0!960.0D0 !0.0D0  
  
  call mpi_comm_rank(mpi_comm_world,myrank,ierror)
  call mpi_comm_size(mpi_comm_world,nproc,ierror)
  
  !master reads input and decomposes the domain  
  if (myrank==master) then
    !read input parameters
    call input_control 
    
    !read geometries
    call input_geom
    
    !intial condition
    call input_ic	
    
    !initialize output time series (?deprecated)
    call init_TS	
    !this should be standarized to catch every boundary condition
    !that may require external information          
  end if
  
  !read friction matrix if neccesary
  call input_friction  
  !decompose the domain, distribute parameters to everyone
  
  call decomp_2d
  
  call printq0
  
  !up to this point
  !everyone should have their one 'global variables'(x,y,z,h,u,v,nbx,nby,etc)
  !so everyone processes their metrics and adimensionalize independently
  !communication has to be made through file BCS.f90, defining a new bc in the vector CB
  !master then decides timestep..
  
!--------------------------------------
  !adimension q and bathy
  call adimension

  !Calcular las metricas
  call metrics
  
  !Crea Coordenadas Xi,Eta
  allocate(coordxi(Nbx),coordeta(Nby))
  coordxi(1)=dxi/2.0D0
  coordeta(1)=deta/2.0D0
  do i=2,Nbx
    coordxi(i)=coordxi(i-1)+dxi
  end do
  do i=2,Nby
    coordeta(i)=coordeta(i-1)+deta
  end do

  !Angulo normales a los bordes con respecto a un eje

  call angulo
!   call VyC(qold_global(3,Nbx,Nby))
  allocate(S1_global(Nbx,Nby),S2_global(Nbx,Nby))
  do i=1,Nbx; do j=1,Nby
	  U1=qold_global(2,i,j)*xi_global(1,i,j)+qold_global(3,i,j)*xi_global(2,i,j)
	  U2=qold_global(2,i,j)*eta_global(1,i,j)+qold_global(3,i,j)*eta_global(2,i,j)
	  S1_global(i,j)=abs(U1)+C_global(i,j)*sqrt(xi_global(1,i,j)**2+xi_global(2,i,j)**2)
	  S2_global(i,j)=abs(U2)+C_global(i,j)*sqrt(eta_global(1,i,j)**2+eta_global(2,i,j)**2)
  end do; end do

  if ( (flagxi0.eq.1).or.(flagxiN.eq.1).or.(flageta0.eq.1).or.(flagetaN.eq.1) )then
    call stability_celerities_boundary_init
  end if

  print*, 'Simulacion Incializada'

END SUBROUTINE init

subroutine printq0
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
  
  write(filename,'("results/grids/grid",I3.3,"_",I3.3,".dat")') coords(1),coords(2)
  open(unit=myrank+100,file=filename)
  write(unit=myrank+100,fmt='( I3.3, "  Nbx" )') Nbx
  write(unit=myrank+100,fmt='( I3.3, "  Nbx" )') Nby
  write(unit=myrank+100,fmt='( I3.3, "  si" )') si
  write(unit=myrank+100,fmt='( I3.3, "  ei" )') ei
  write(unit=myrank+100,fmt='( I3.3, "  sj" )') sj
  write(unit=myrank+100,fmt='( I3.3, "  ej" )') ej
  write(unit=myrank+100,fmt='( I3.3, "  coord(1)" )') coords(1)
  write(unit=myrank+100,fmt='( I3.3, "  coord(2)" )') coords(1)
  write(unit=myrank+100,fmt='( I3.3, "  dims(1)" )') dims(1)
  write(unit=myrank+100,fmt='( I3.3, "  dims(2)" )') dims(1)
  do i=1,Nbx
    do j=1,Nbx
      write(unit=myrank+100,fmt='(6(E15.5,2X))') x_global(i,j),y_global(i,j), z_global(i,j)
    end do
  end do
  close(unit=myrank+100)
  
  !copy gridX,gridY and gridZ
  if (myrank==0) then
    do i=1,ngrids
      if( (batiopt(i)==0).or.(batiopt(i)==1) )then
	call system('cp '//batiname(i,1)//' results/.')
	call system('cp '//batiname(i,2)//' results/.')
	call system('cp '//batiname(i,3)//' results/.')
      else if( (batiopt(i)==2).or.(batiopt(i)==3) )then
	call system('cp '//batiname(i,1)//' results/.')
      end if
    end do
  end if
end subroutine printq0

SUBROUTINE ADIMENSION
  !Function that aplies the adimensionalization to the initial conditions
  !Funcion que aplica la adimensionalizacion a las condiciones inciales y a todo
  USE global_variables
  USE geometries
  implicit none
  integer :: i,j
  allocate(V_global(Nbx,Nby),C_global(Nbx,Nby),VC(Nbx,Nby))
  ! real (kind=8)::U1,U2
  do i=1,Nbx; do j=1,Nby
    x_global(i,j)=x_global(i,j)/L
    y_global(i,j)=y_global(i,j)/L
    z_global(i,j)=z_global(i,j)/H
    qold_global(1,i,j)=qold_global(1,i,j)/H
    qold_global(2,i,j)=qold_global(2,i,j)/U
    qold_global(3,i,j)=qold_global(3,i,j)/U
    V_global(i,j)=sqrt((qold_global(2,i,j))**2.0D0+(qold_global(3,i,j))**2.0D0)
    C_global(i,j)=sqrt(qold_global(1,i,j)/FR2)
    VC(i,j)=V_global(i,j)+C_global(i,j)			
  end do; end do
END SUBROUTINE ADIMENSION
