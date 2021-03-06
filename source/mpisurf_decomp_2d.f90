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

  
  !mpi knows how to split comm_world better into a cart topo
  call mpi_dims_create(nproc,ndim,dims,ierror)
  call mpi_cart_create(mpi_comm_world,ndim,dims,isperiodic,&
    reorder,comm2d,ierror)
  call mpi_comm_rank(comm2d,myrank2d,ierror)
  !get the coords for this core
  call mpi_cart_get(comm2d,ndim,dims,isperiodic,coords,ierror)
  coords2=coords
  
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
  !bcast is simpler (way simpler) to (human) read!
  
  !bcast input control parameters
  
  call bcast_input_control
  
  !distribute friction parameters
  call bcast_friction
  
  !distribuite geometries and initial condition
  call bcast_initq_grids

  !properly assign boundary conditions...
  !in a separate file (a lot of lines)
  call bcast_bcs
! print*,'asdf',myrank      
  !bcast gauge points
!   print*,'a',myrank
  call bcast_gauges
!   print*,'b',myrank
end subroutine decomp_2d

subroutine bcast_input_control
  !routine to broadcast input_control.f90 parameters 
  !it is good to everyone know about this
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
  Nbx=ei-si+1
  Nby=ej-sj+1
  
  !allocate/initialize everything,for everyone
  !in the same order as input_control.f90
  !careful with datatypes!  
  call mpi_bcast(caso,	1,	mpi_integer,master,comm2d,ierror)
  call mpi_bcast(tinit,1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(treal,1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(t ,	1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(dt,	1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(dtreal,1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(tfinal,1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(cfl,	1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(dxi,	1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(deta,	1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(L,	1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(H,	1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(U,	1,	mpi_double_precision,master,comm2d,ierror)  
  call mpi_bcast(CB_real,4,	mpi_integer,master,comm2d,ierror)
  call mpi_bcast(dit,	1,	mpi_integer,master,comm2d,ierror)
  if (dit==-1) then
    call mpi_bcast(dtout,1,	mpi_double_precision,master,comm2d,ierror)
  end if
  call mpi_bcast(kappa,1,	mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(rk,	1,	mpi_integer,master,comm2d,ierror)
  call mpi_bcast(mmopt,1,	mpi_integer,master,comm2d,ierror)
  call mpi_bcast(fopt,	1,	mpi_integer,master,comm2d,ierror) 
  call mpi_bcast(outopt,1,	mpi_integer,master,comm2d,ierror)
  call mpi_bcast(outdir,255,	mpi_character,master,comm2d,ierror)
  call mpi_bcast(indir,255,	mpi_character,master,comm2d,ierror)
  call mpi_bcast(Fr2,	1,	mpi_double_precision,master,comm2d,ierror)  
end subroutine

subroutine bcast_friction
  use global_variables
  use mpi
  use mpi_surf
  use multigrid_surf
  implicit none
  integer::level
  if (fopt==0) then
      Cf=0    !no need to bcast, each proc assigns a 0
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
     
  x_global=geom(myzone)%X(si:ei,sj:ej)  
  y_global=geom(myzone)%Y(si:ei,sj:ej)
  z_global=geom(myzone)%Z(si:ei,sj:ej)  
  qold_global(1,:,:)=initq(myzone)%H(si:ei,sj:ej)
  qold_global(2,:,:)=initq(myzone)%U(si:ei,sj:ej)
  qold_global(3,:,:)=initq(myzone)%V(si:ei,sj:ej)

  deallocate(geom,initq)
end subroutine bcast_initq_grids

subroutine bcast_gauges
  use TimeSeries
  use geometries, only: x_global, y_global
  use global_variables, only: outdir
  use mpi_surf
  use multigrid_surf
  use mpi
  
  implicit none
  integer :: i,aloc
  character(len=100)::intchar
  !idea: use m1_temp to see if the point lies in mygrid

  call mpi_bcast(Nts,1, mpi_integer,master,comm2d,ierror)
  
  !let other procs allocate ids,x0,y0, and m1_temp
  if (myrank/=master) then
    allocate(x0(Nts),y0(Nts),id0(Nts), m1_temp(2,Nts))
  end if

  !bcast this thing  
  call mpi_bcast(m1_temp,2*Nts, mpi_integer,master,comm2d,ierror)
  call mpi_bcast(x0,Nts, mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(y0,Nts, mpi_double_precision,master,comm2d,ierror)
  call mpi_bcast(id0,Nts, mpi_integer,master,comm2d,ierror)  
  
  call mpi_barrier(comm2d,ierror)
  !allocate position indices and gaugeflags
  allocate(m1(2,Nts),gaugeflag(Nts))
  !now find the indices

  write(intchar,'("P",I3.3,"R",I3.3)') nproc,myrank
  open(unit=10+myrank,file=trim(outdir)//'/timeseries/gaugeflag'//trim(adjustl(intchar))//'.txt')
  DO i=1,Nts    
    !im assuming rectangular grid
    !otherwise the criteria is just to check if 
    !the points is inside the rectangle on which
    !the grid is inscribed        
    if (m1_temp(1,i)>=si .and. m1_temp(1,i)<=ei .and. &
	m1_temp(2,i)>=sj .and. m1_temp(2,i)<=ej) then	
	gaugeflag(i) = .true.
	m1(:,i) = minloc((x_global-x0(i))*(x_global-x0(i))+(y_global-y0(i))*(y_global-y0(i)))    
    else
	gaugeflag(i) = .false.
    end if    
      write(10+myrank,*) gaugeflag(i)
 END DO
  
  close(10+myrank)
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
end subroutine decomp1d
