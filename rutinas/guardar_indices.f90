SUBROUTINE guardar_indices
!rutina para guardar los indices de las celdas del dominio que contienen 
!genera un archivo con las columnas
!id,x,y,ix,iy,d
!id=numero de sonda
!x,y las coordenadas
!ix,iy indices de las coordenadas
!d es la distancia entre el punto exacto y el centro de la celda
USE global_variables
USE geometries
USE LagrangianPT
USE TimeSeries
implicit none

real (kind=8) :: dr !la distancia entre x0,y0 y la celda (x,y)(m1)
integer	:: i,j, ent, itent, ent1
real (kind=8)	:: ind, dec, time, dit1
real (kind=8)	:: ind1, dec1
logical	:: lexist, lexistT
 character(len=1100):: filename, filenameT,fgauges
 character(len=10)::number, numbercaso
 character(len=1000)::intchar, ncaso
 character(len=1100):: filename_LPTx, filename_LPTy, filename_LPTu, filename_LPTv
!-------------------------
logical	:: lexistPh, lexistPu, lexistPv, lexistT_eta, lexist_LPTx, lexist_LPTy, lexist_LPTu, lexist_LPTv 
 character(len=1100):: filenamePh,filenamePu,filenamePv,filenameT_eta, path
logical	:: lexist_param
 character(len=1100):: filename_param
logical	:: lexist_TSxy
 character(len=1100):: filename_TSxy
logical	:: lexist_H, lexist_U, lexist_V
 character(len=1100):: filename_H, filename_U, filename_V
 real (kind=8) :: h_g,u_g,v_g,z_g

path='results/'
fgauges='data/indices.dat'


do i=1,Nts
		!1.-Encontrar indices del vertice superior derecho de la celda que contiene al punto i
! 		m1=minloc((x_global-x0(i))*(x_global-x0(i))+(y_global-y0(i))*(y_global-y0(i)), &
! 			  mask=(x_global-x0(i)).ge.0 .and. &
! 			       (y_global-y0(i)).ge.0 )	
		!nope,mejor encontrar LA celda más cercana, pues son centradas y representan con ese valor a toda esa region
		m1=minloc((x_global-x0(i))*(x_global-x0(i))+(y_global-y0(i))*(y_global-y0(i)))
	    
		inquire(FILE=fgauges, EXIST=lexist)
		print*,lexist
		if (.NOT. lexist) then
		  open(20,file=fgauges,status='new',action='write')
! 		  pause
		else    
		  if (it==0.0D0 .and. i==1) then
		    open(20,file=fgauges,status='replace',action='write')
		  else
		    open(20,file=fgauges,status='old',action='write', position='append')
		  end if
		end if
		!se podria guardar id,x,y, en otro archivo aparte
		dr=sqrt((x_global(m1(1),m1(2))-x0(i))*(x_global(m1(1),m1(2))-x0(i))+(y_global(m1(1),m1(2))-y0(i))*(y_global(m1(1),m1(2))-y0(i)))
		write(20,21) id0(i),x0(i),y0(i),int(m1(1)),int(m1(2)),dr
		21 format(I, 1X, F15.5, 1X, F15.5, 1X, I, 1X, I, 1X, F15.5 /)
		close(20)
		
end do
END SUBROUTINE guardar_indices


