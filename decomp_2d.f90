!subroutines in this file:
!	decomp_2d: 	is like the 'mainfile' for decomposing 
!	bcast_input_control
!	bcast_friction
!	decomp1d: splits range(N) into nprocs in a 'balanced' way
!	set_bcs: assign boundary conditions...distribuite input time series and give
!		a key number to communication boundaries

subroutine decomp_2d
  use global_variables
  use geometries
  use mpi
  use mpi_surf
  use multigrid_surf
  implicit none
  
  !split comm_world by zones (not yet)
  
  
  !split each zone in a topology  
  !catch periodic boundaries for a better topology
  !get dimensions of topology
  
  
  isperiodic(1)=.false.
  isperiodic(2)=.false.
  reorder=.true.  
  if (CB_real(1)==2.and.CB_real(2)==2) then
    isperiodic(1)=.true.
  else if (CB_real(1)==2.and.CB_real(2)/=2) then
    print*,'wrong boundary in bc(1) or bc(2) (periodic)'
    call mpi_abort(comm2d,100,ierror)
  else if (CB_real(1)/=2.and.CB_real(2)==2) then
    print*,'wrong boundary in bc(1) or bc(2) (periodic)'
    call mpi_abort(comm2d,100,ierror)
  end if
  
  if (CB_real(3)==2.and.CB_real(4)==2) then
    isperiodic(2)=.true.
  else if (CB_real(3)==2.and.CB_real(4)/=2) then
    print*,'wrong boundary in bc(3) or bc(4) (periodic)'
    call mpi_abort(comm2d,100,ierror)
  else if (CB_real(3)/=2.and.CB_real(4)==2) then
    print*,'wrong boundary in bc(3) or bc(4) (periodic)'
    call mpi_abort(comm2d,100,ierror)
  end if

  
  !mpi knows how to split comm_world better into a topo
  call mpi_dims_create(nproc,ndim,dims,ierror)
  call mpi_cart_create(mpi_comm_world,ndim,dims,isperiodic,&
    reorder,comm2d,ierror)
  call mpi_comm_rank(comm2d,myrank2d,ierror)
  !get the coords for this core
  call mpi_cart_get(comm2d,ndim,dims,isperiodic,coords,ierror)
  
  !get my neighbors
  !left/right & back/front as when looking to the matrix
  !increasing indices to the upper right corner 
  !ex:  1,0	1,1
  !	0,0	0,1   , left(1,1)=1,0;  back(1,1)=0,1
  call mpi_cart_shift(comm2d,0,shift,myback,myfront,ierror) 
  call mpi_cart_shift(comm2d,1,shift,myleft,myright,ierror) 
  
  
  
  !now, broadcast  everything...or almost..
  !(still think sendrcv would be better, at least in terms of memory
  !but i checked experimentally, and it was worse, i may have done it wrong)
  !bcast is simpler (way simpler) to read!
  
  !bcast input control parameters
  call bcast_input_control
  
  !distribute friction parameters
  call bcast_friction
  
  !distribuite geometries and initial condition
  call bcast_initq_grids
  
  !properly assign boundary conditions...
  call set_bcs
end subroutine decomp_2d

subroutine bcast_input_control
  use mpi
  use mpi_surf
  use multigrid_surf
  use global_variables
  implicit none  
  
  !once i get the overlapping thing, these two bcast should be done
  !in the zone-splitted commm
  call mpi_bcast(ngrids,1,mpi_integer,master,comm2d,ierror)
  if (myrank/=master) then
    allocate(nxi(ngrids),neta(ngrids))
  end if
  
  call mpi_bcast(nxi,ngrids,mpi_integer,master,comm2d,ierror)
  call mpi_bcast(neta,ngrids,mpi_integer,master,comm2d,ierror)
  call decomp1d(nxi(myzone),dims(1),coords(1),si,ei)
  call decomp1d(neta(myzone),dims(2),coords(2),sj,ej)    
  !allocate/initialize everything,for everyone
  !in the same order as input_control.f90
  call mpi_bcast(caso,1,mpi_integer,master,comm2d,ierror)
  call mpi_bcast(tinit,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(tfinal,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(cfl,1,mpi_double_precision,master,comm2d,ierror)
  Nbx=ei-si+1
  Nby=ej-sj+1
  call mpi_bcast(dxi,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(deta,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(L,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(H,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(U,1,mpi_double_precision,master,comm2d,ierror)  
  call mpi_bcast(CB_real,4,mpi_integer,master,comm2d,ierror)
  call mpi_bcast(dit,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(kappa,1,mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(rk,1,mpi_integer,master,comm2d,ierror)
  call mpi_bcast(mmopt,1,mpi_integer,master,comm2d,ierror)
  call mpi_bcast(fopt,1,mpi_integer,master,comm2d,ierror) 
  call mpi_bcast(outopt,1,mpi_integer,master,comm2d,ierror)
  call mpi_bcast(Fr2,1,mpi_double_precision,master,comm2d,ierror)  
end subroutine

subroutine bcast_friction
  use global_variables
  use mpi
  use mpi_surf
  use multigrid_surf
  implicit none
  integer::level
  
  if (fopt==0) then
      Cf=0    !dont need to bcast, each proc assigns a 0
      Coef=0.0D0
  else if (fopt==1) then
    call mpi_bcast(fM,1,mpi_integer,master,comm2d,ierror)
    call mpi_bcast(Cf,1,mpi_integer,master,comm2d,ierror)
    if (fM==1) then
      call mpi_bcast(Coef,1,mpi_double_precision,comm2d,ierror)
    else if (fM==2) then      
      if (myrank/=master) then !allocate buffer for non_master
	allocate(buffMcoef(ngrids))	
      end if
	!i should learn how to bcast a structure
	do level=1,ngrids
	  if (myrank/=master) then
	    allocate(buffMcoef(level)%bMCoef(nxi(level),neta(level)))
	  end if
	  call mpi_bcast(buffMcoef(level)%bMcoef(nxi(level),neta(level)),&
	    nxi(level)*neta(level),mpi_double_precision,master,comm2d,ierror)	  
	end do	
      !assign Mcoef respectively
      allocate(Mcoef(Nbx,Nby))
      MCoef(:,:)=buffMcoef(myzone)%bMcoef(si:ei,sj:ej)
      !destroy the buffer
      deallocate(buffMcoef)     
    end if
  end if  
end subroutine

subroutine bcast_initq_grids
  use global_variables
  use geometries
  use mpi
  use mpi_surf
  use multigrid_surf
  implicit none
  integer::level,i,j
  if (myrank/=master) then
    allocate(geom(ngrids),initq(ngrids))
  end if
  allocate(x_global(Nbx,Nby),y_global(Nbx,Nby),z_global(Nbx,Nby))
  allocate(qold_global(3,Nbx,Nby))
  do level=1,ngrids
    if (myrank/=master) then
      allocate(geom(level)%X(nxi(level),neta(level)), &
	geom(level)%Y(nxi(level),neta(level)), &
	geom(level)%Z(nxi(level),neta(level)), &
	initq(level)%H(nxi(level),neta(level)), &
	initq(level)%U(nxi(level),neta(level)), &
	initq(level)%V(nxi(level),neta(level)))    
    end if
    !bcast geom       
    call mpi_bcast(geom(level)%X,nxi(level)*neta(level), &
      mpi_double_precision,master,comm2d,ierror)
    call mpi_bcast(geom(level)%Y,nxi(level)*neta(level), &
      mpi_double_precision,master,comm2d,ierror)
    call mpi_bcast(geom(level)%Z,nxi(level)*neta(level), &
      mpi_double_precision,master,comm2d,ierror)
    
    call mpi_bcast(initq(level)%H,nxi(level)*neta(level), &
      mpi_double_precision,master,comm2d,ierror)
    call mpi_bcast(initq(level)%U,nxi(level)*neta(level), &
      mpi_double_precision,master,comm2d,ierror)      
    call mpi_bcast(initq(level)%V,nxi(level)*neta(level), &
      mpi_double_precision,master,comm2d,ierror)    
  end do  
     
  x_global(:,:)=geom(myzone)%X(si:ei,sj:ej)  
  y_global(:,:)=geom(myzone)%Y(si:ei,sj:ej)
  z_global(:,:)=geom(myzone)%Z(si:ei,sj:ej)  
  qold_global(1,:,:)=initq(myzone)%H(si:ei,sj:ej)
  qold_global(2,:,:)=initq(myzone)%U(si:ei,sj:ej)
  qold_global(3,:,:)=initq(myzone)%V(si:ei,sj:ej)

  deallocate(geom,initq)
end subroutine bcast_initq_grids

subroutine set_bcs
  !set CB vector and scatter input time series if neccesary
  use global_variables
  use mpi_surf
  use multigrid_surf
  use mpi
  
  CB(1)=-1
  CB(2)=-1
  CB(3)=-1
  CB(4)=-1  
  
  if (myback==mpi_proc_null) then 
    CB(1)=CB_real(1)
  end if
  
  if (myfront==mpi_proc_null) then
    CB(2)=CB_real(2)
  end if
  
  if (myleft==mpi_proc_null) then
    CB(3)=CB_real(3)
  end if
  
  if (myright==mpi_proc_null) then
    CB(4)=CB_real(4)
  end if
  
  !catch periodic boundaries
  if (CB(1)==2) then
    CB(1)=-1
    CB(2)=-1
  end if
  
  if (CB(3)==2) then
    CB(3)=-1
    CB(4)=-1    
  end if
  
  
end subroutine
subroutine decomp1d(n,p,rank,s,e) !based in the book of Gropp,1999
  !recibe el numero total n y el numero de bloques p
  !rank es el indice del que queremos saber
  !los numeros s y e
  !s es el indice inferior
  !e es el indice superior
  !la forma mas facil es ver que pasa con s_rank, y e=s_{rank+1}-1
  !para los diferentes casos: rank<deficit o rank>=deficit
  !si rank<deficit s=s+r, sino s=s+deficit....esa es la diferencia
  !la idea es pasar el deficit a los primeros procesadores y asi balancear la cosa
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
