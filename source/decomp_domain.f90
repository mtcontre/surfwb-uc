subroutine decomp_domain
  use mpi
  use mpi_surf
  use global_variables
  use geometries
  implicit none
  !get the starting indexes in the i-j direction (si,ei, sj,ej)
  !to split each matrix into pieces
  
  call decomp1d(nxi, dims(1), coords(1), si, ei)
  call decomp1d(neta, dims(2), coords(2), sj, ej)
  
  Nbx = ei-si+1
  Nby = ej-sj+1
  
  !grids+bati
  allocate(x_global(nbx,nby), y_global(nbx,nby), z_global(nbx,nby))  
  x_global(:,:) = x_buff(si:ei,sj:ej)
  y_global(:,:) = y_buff(si:ei,sj:ej)
  z_global(:,:) = z_buff(si:ei,sj:ej)  
  deallocate(x_buff, y_buff, z_buff)
  
  !initq
  allocate(qold_global(3,nbx,nby),qnew_global(3,Nbx,Nby),qreal_global(3,Nbx,Nby))
  qold_global(:,:,:) = q_buff(:, si:ei, sj:ej)
  !initialize qout for 1st step					!antes estaba en outputmat..pero se perdia
  qnew_global(:,:,:) = qold_global(:,:,:)
  deallocate(q_buff)
  
  !friction
  allocate(Mcoef(nbx,nby))
  Mcoef(:,:) = fric_buff(si:ei, sj:ej)
  deallocate(fric_buff)
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
end subroutine decomp1d