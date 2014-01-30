subroutine couplingbc_xi0
  !load coupling boundary conditions for xi=0
  !...for now only reads using 
  use global_variables
  use couplingbc
  use multigrid_surf

  implicit none
  character(len=3) :: sep
  integer:: nvar,ivar,ix,indt !nvar={1:z,2:h,3:u,4:v}
  nvar=4!number of variables to read: (z,h,u,v)
  print*,'Reading coupling BC for xi_0'
  !1: first group (row,column,whatever) of ghost cells
  open(unit=99,file=fnamexi0(1))
  if (optxi0g1==0) then
    print*,'xi0g1 boundary interpolation: piecewise constant'
  else
    print*,'xi0g1 boundary interpolation: linear'
  end if 
  allocate(bufqxi0g1(nvar,neta(1),nt_xi0g1+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(99,*)((bufqxi0g1(ivar,ix,indt),indt=1,nt_xi0g1+1),ix=1,neta(1))
    if (ivar<=nvar-1) then
      read(99,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=99)
  !2: second group of ghost cells
  open(unit=100,file=fnamexi0(2))
  if (optxi0g2==0) then
    print*,'xi0g2 boundary interpolation: piecewise constant'
  else
    print*,'xi0g2 boundary interpolation: linear'
  end if
  allocate(bufqxi0g2(nvar,neta(1),nt_xi0g2+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(100,*)((bufqxi0g2(ivar,ix,indt),indt=1,nt_xi0g2+1),ix=1,neta(1))
    if (ivar<=nvar-1) then
      read(100,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=100)
end subroutine couplingbc_xi0

subroutine couplingbc_xiN
  !load coupling boundary conditions for xi=0
  !2 
  use global_variables
  use couplingbc
  use multigrid_surf
  implicit none
  character(len=3) :: sep
  integer:: nvar,ivar,ix,indt !nvar={1:z,2:h,3:u,4:v}
  nvar=4!number of variables to read: (z,h,u,v)
  print*,'Reading coupling BC for xi_N'
  !1: first group (row,column,whatever) of ghost cells
  open(unit=99,file=fnamexiN(1))

  if (optxiNg1==0) then
    print*,'xiNg1 boundary interpolation: piecewise constant'
  else
    print*,'xiNg1 boundary interpolation: linear'
  end if  
  allocate(bufqxiNg1(nvar,neta(1),nt_xiNg1+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(99,*)((bufqxiNg1(ivar,ix,indt),indt=1,nt_xiNg1+1),ix=1,neta(1))
    if (ivar<=nvar-1) then
      read(99,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=99)
  !2: second group of ghost cells
  open(unit=100,file=fnamexiN(2))

  if (optxiNg2==0) then
    print*,'xiNg2 boundary interpolation: piecewise constant'
  else
    print*,'xiNg2 boundary interpolation: linear'
  end if
  allocate(bufqxiNg2(nvar,neta(1),nt_xiNg2+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(100,*)((bufqxiNg2(ivar,ix,indt),indt=1,nt_xiNg2+1),ix=1,neta(1))
    if (ivar<=nvar-1) then
      read(100,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=100)
end subroutine couplingbc_xiN

subroutine couplingbc_eta0
  !load coupling boundary conditions for xi=0
  !2 
  use global_variables
  use couplingbc
  use multigrid_surf
  implicit none
  character(len=3) :: sep
  integer:: nvar,ivar,ix,indt !nvar={1:z,2:h,3:u,4:v}
  nvar=4!number of variables to read: (z,h,u,v)
  print*,'Reading coupling BC for eta0'
  !1: first group (row,column,whatever) of ghost cells
  open(unit=99,file=fnameeta0(1))
  if (opteta0g1==1) then
    print*,'eta0g1 boundary interpolation: piecewise constant'
  else
    print*,'eta0g1 boundary interpolation: linear'
  end if
  allocate(bufqeta0g1(nvar,nxi(1),nt_eta0g1+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(99,*)((bufqeta0g1(ivar,ix,indt),indt=1,nt_eta0g1+1),ix=1,nxi(1))
    if (ivar<=nvar-1) then
      read(99,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=99)
  !2: second group of ghost cells
  open(unit=100,file=fnameeta0(2))
  if (opteta0g2==0) then
    print*,'eta0g2 boundary interpolation: piecewise constant'
  else
    print*,'eta0g2 boundary interpolation: linear'
  end if
  allocate(bufqeta0g2(nvar,nxi(1),nt_eta0g2+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(100,*)((bufqeta0g2(ivar,ix,indt),indt=1,nt_eta0g2+1),ix=1,nxi(1))
    if (ivar<=nvar-1) then
      read(100,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=100)  
end subroutine couplingbc_eta0

subroutine couplingbc_etaN
  !load coupling boundary conditions for xi=0
  !2 
  use global_variables
  use couplingbc
  use multigrid_surf
    use mpi_surf
  implicit none
  character(len=3) :: sep
  integer:: nvar,ivar,ix,indt !nvar={1:z,2:h,3:u,4:v}
  nvar=4!number of variables to read: (z,h,u,v)
  print*,'Reading coupling BC for etaN'  
  !1: first group (row,column,whatever) of ghost cells
  open(unit=99,file=fnameetaN(1))

  if (optetaNg1==0) then
    print*,'etaNg1 boundary interpolation: piecewise constant'
  else
    print*,'etaNg1 boundary interpolation: linear'
  end if
  allocate(bufqetaNg1(nvar,nxi(1),nt_etaNg1+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(99,*)((bufqetaNg1(ivar,ix,indt),indt=1,nt_etaNg1+1),ix=1,nxi(1))
    if (ivar<=nvar-1) then
      read(99,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=99)
  
  !2: second group of ghost cells
  open(unit=100,file=fnameetaN(2))

  if (optetaNg2==0) then
    print*,'etaNg2 boundary interpolation: piecewise constant'
  else
    print*,'etaNg2 boundary interpolation: linear'
  end if
  allocate(bufqetaNg2(nvar,nxi(1),nt_etaNg2+1))  !(nvars,ny,nt+1)
  do ivar=1,nvar
    read(100,*)((bufqetaNg2(ivar,ix,indt),indt=1,nt_etaNg2+1),ix=1,nxi(1))
    if (ivar<=nvar-1) then
      read(100,*) !this line is just to separate groups of variables
      !actually usefull when reviewing data 'by eye'
    end if
  end do
  close(unit=100)
  
  
end subroutine couplingbc_etaN

! CONDICION DE BORDE PARA EL PRIMER GRUPO (COLUMNA O FILA) DE CELDAS FANTASMA
! Codigo por Jose Galaz Mora (jdgalaz@uc.cl),2013
! para el Depto. de Ing. Hidr. y Ambiental de la Esc. de Ing. de la Pontificia Universidad Catolica 
! 
! 1º Formato de input: Matriz q(t) definida por encabezado+matrizq1+matrizq2+matrizq3
! 
! 100	nt    numero de pasos de tiempo sin contar el inicial
! 0.1	dt    intervalo de tiempo
! 1	opt   dice si utilizar escalones (=1) o interpolación lineal (otro caso)
! ***	sep   separador de bloque,puede ser cualquier cosa=caracter que separa q(1) de q(2) y q(2) de q(3)	      
! q(1,y=y0,t=0)		q(1,y=y0,t=dt)		q(1,y=y0,t=2dt)		q(1,y=y0,nt*dt)
! q(1,y=y1,t=0)		q(1,y=y1,t=dt)		q(1,y=y1,t=2dt)		q(1,y=y1,nt*dt) --
! q(1,y=y2,t=0)		q(1,y=y2,t=dt)		q(1,y=y2,t=2dt)		q(1,y=y2,nt*dt) --
! ...
! q(1,y=y(neta),t=0)	q(1,y=y(neta),t=dt)	q(1,y=y(neta),t=2dt)	q(1,y=y(neta),nt*dt)
! ***
! q(2,y=y0,t=0)		q(2,y=y0,t=dt)		q(2,y=y0,t=2dt)		q(2,y=y0,nt*dt)
! q(2,y=y1,t=0)		q(2,y=y1,t=dt)		q(2,y=y1,t=2dt)		q(2,y=y1,nt*dt) --
! q(2,y=y2,t=0)		q(2,y=y2,t=dt)		q(2,y=y2,t=2dt)		q(2,y=y2,nt*dt) --
! ...
! q(2,y=y(neta),t=0)	q(2,y=y(neta),t=dt)	q(2,y=y(neta),t=2dt)	q(2,y=y(neta),nt*dt)
! ***
! q(3,y=y0,t=0)		q(3,y=y0,t=dt)		q(3,y=y0,t=2dt)		q(3,y=y0,nt*dt)
! q(3,y=y1,t=0)		q(3,y=y1,t=dt)		q(3,y=y1,t=2dt)		q(3,y=y1,nt*dt) --
! q(3,y=y2,t=0)		q(3,y=y2,t=dt)		q(3,y=y2,t=2dt)		q(3,y=y2,nt*dt) --
! ...
! q(3,y=y(neta),t=0)	q(3,y=y(neta),t=dt)	q(3,y=y(neta),t=2dt)	q(3,y=y(neta),nt*dt)
! 
! 2º Importar las variables
!   Un entero con el valor de nt
!   Un real con el valor de dt
!   Matriz de rango 3: qg1(i,j,k): q-ghost-1(variable,celda,tiempo)
! 		      i={1,2,3,4}:(z,h,u,v)-->variable
! 		      j=1,2,3,..neta
! 		      k={1,2,3,..,nt}
!   por lo tanto: 
!   
!   en variables globales: definir nt_xi0g1, optxi0g1, dt_xi0g1 ,qxi0g1(:,:,:) (xi0g1 por xi=0 ghost 1)
!   
!   (no defino el vector tiempo, ya que con nt y dt tengo toda la información necesaria)
!   
! 2ºAsignar para un tiempo t* que no necesariamente pertenece a t!!
!    ---interpolar linealmente? (enfoque continuo)---> por el runge kutta??
!    ---aproximar por el menor valor más cercano??---> enfoque vol.finitos...MEJOR este, sino suavizo las discontinuidades
!    queda entonces 'idea': q(t)=f(int(t)), por escalones
! 3ºPasar estas variables al qt

!!!!!!!!!!!bcxi01 y bcxi02 son para las celdas -2 y -1
!!!!!!!!!! bcxiN1 y bcxiN2 son para las celdas Nby+1 y Nby+2
!!!!!!!!!!es decir, 'siguen el orden en que se leen'