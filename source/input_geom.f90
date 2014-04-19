SUBROUTINE input_geom
  use global_variables
  use geometries
  use multigrid_surf
  use mpi_surf
  implicit none

  real (kind=8), dimension(:), allocatable :: xaux, yaux, zaux
  integer	:: i,j,level !level=grid level
  allocate(geom(ngrids))
  !initialize each grid
  !maybe, when using overlapped grids
  !i should split comm_world in ngrids communicators
  !so each master reads each input file (which are the heaviest)
  do level=1,ngrids
    allocate(geom(level)%X(nxi(level),neta(level)), &
      geom(level)%Y(nxi(level),neta(level)), &
      geom(level)%Z(nxi(level),neta(level)))      
    SELECT CASE (int(batiopt(level)))
      CASE(0)
	open(unit=2,file=batiname(level,1),form='unformatted')
	read(2) ((geom(level)%X(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
	
	open(unit=2,file=batiname(level,2),form='unformatted')
	read(2) ((geom(level)%Y(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)	
	
	open(unit=2,file=batiname(level,3),form='unformatted')
	read(2) ((geom(level)%Z(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
      CASE(1)
	open(unit=2,file=batiname(level,1))
	read(2,*) ((geom(level)%X(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
	
	open(unit=2,file=batiname(level,2))
	read(2,*) ((geom(level)%Y(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
	
	open(unit=2,file=batiname(level,3))
	read(2,*) ((geom(level)%Z(i,j),j=1,neta(level)),i=1,nxi(level))
	close(unit=2)
      CASE(2)	
	allocate(xaux(nxi(level)*neta(level)),yaux(nxi(level)*neta(level)), zaux(nxi(level)*neta(level)))
	open(unit=2,file=batiname(level,1))
	Do i=1,nxi(level)*neta(level)
	  read(2,*) xaux(i), yaux(i), zaux(i)
	End Do
	close(unit=2)
	do i=1,nxi(level),1; do j=1,neta(level),1
	  geom(level)%X(i,j)=xaux(j+(i-1)*neta(level))
	  geom(level)%Y(i,j)=yaux(j+(i-1)*neta(level))
	  geom(level)%Z(i,j)=zaux(j+(i-1)*neta(level))
	end do; end do
	deallocate(xaux,yaux,zaux)
      CASE(3)
	allocate(xaux(nxi(level)*neta(level)),yaux(nxi(level)*neta(level)), zaux(nxi(level)*neta(level)))
	open(unit=2,file=batiname(level,1))
	Do i=1,nxi(level)*neta(level)
	  read(2,*) xaux(i), yaux(i), zaux(i)
	End Do
	close(unit=2)
	do j=1,neta(level),1; do i=1,nxi(level),1
	  geom(level)%X(i,j)=xaux(i+(j-1)*nxi(level))
	  geom(level)%Y(i,j)=yaux(i+(j-1)*nxi(level))
	  geom(level)%Z(i,j)=zaux(i+(j-1)*nxi(level))
	end do; end do
	deallocate(xaux,yaux,zaux)
    END SELECT
  end do
  
  !for debugging..show first 3row/columns
!   print*,'X------------'
!   level=1
!   do i =1,3
!     print*,geom(level)%X(i,1),geom(level)%X(i,2),geom(level)%X(i,3)
!   end do
! !   pause
!   
!   print*,'Y------------'
!   do i =1,3
!     print*,geom(level)%Y(i,1),geom(level)%Y(i,2),geom(level)%Y(i,3)
!   end do
! !   pause
!   
!   print*,'Z------------'
!   do i =1,3
!     print*,geom(level)%Z(i,1),geom(level)%Z(i,2),geom(level)%Z(i,3)
!   end do
END SUBROUTINE input_geom
