SUBROUTINE outputgauges(time)

USE global_variables
USE geometries
USE LagrangianPT
USE TimeSeries
implicit none

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

!! Coordenadas Puntos
! x_1,x_2,y_1,y_2,,f11,f12,fs21,f22
!Dimensionalization of the variables
do i=1,Nbx; do j=1,Nby
		qreal_global(1,i,j)=qnew_global(1,i,j)*H
		qreal_global(2,i,j)=qnew_global(2,i,j)*U
		qreal_global(3,i,j)=qnew_global(3,i,j)*U
end do; end do
if (time==0.0D0) then
qnew_global=qold_global

end if
write(ncaso,*) caso
ncaso=adjustl(ncaso)
numbercaso=ncaso  
path='results/'
write(ncaso,*) caso
ncaso=(adjustl(ncaso))
fgauges='results/timeseries/gauges.' // trim(ncaso)// '.dat'


do i=1,Nts
		!1.-Encontrar indices del vertice superior derecho de la celda que contiene al punto i
		m1=minloc((x_global-x0(i))*(x_global-x0(i))+(y_global-y0(i))*(y_global-y0(i)), &
			  mask=(x_global-x0(i)).ge.0 .and. &
			       (y_global-y0(i)).ge.0 )	
		
		!3.- Generar los datos interpolados (interpolacion bilineal)
		!------------------h-------------------------------
		call interpxy(m1,x_global,y_global,qreal_global(1,:,:),x0(i),y0(i),h_g)
		call interpxy(m1,x_global,y_global,qreal_global(2,:,:),x0(i),y0(i),u_g)
		call interpxy(m1,x_global,y_global,qreal_global(3,:,:),x0(i),y0(i),v_g)
		call interpxy(m1,x_global,y_global,z_global(:,:),x0(i),y0(i),z_g)
! 		pause
! 
! 		!------------------u-------------------------------
! 		!Obtener la f de las coordenadas de los vertices
! 		f11=qreal_global(2,m1(1)-1,m1(2)-1)
! 		f12=qreal_global(2,m1(1)-1,m1(2))
! 		f21=qreal_global(2,m1(1)	 ,m1(2)-1)
! 		f22=qreal_global(2,m1(1)	 ,m1(2)	)
! 		!Interpolante bilineal
! 		u_g=( f11*(x_2-x0(i))*(y_2-y0(i))+ &
! 		      f21*(x0(i)-x_1)*(y_2-y0(i))+ &
! 		      f12*(x_2-x0(i))*(y0(i)-y_1) + &
! 		      f22*(x0(i)-x_1)*(y0(i)-y_1) )/((x_2-x_1)*(y_2-y_1)) 
! 
! 		!------------------v-------------------------------
! 		!Obtener la f de las coordenadas de los vertices
! 		f11=qreal_global(3,m1(1)-1,m1(2)-1)
! 		f12=qreal_global(3,m1(1)-1,m1(2))
! 		f21=qreal_global(3,m1(1)	 ,m1(2)-1)
! 		f22=qreal_global(3,m1(1)	 ,m1(2)	)
! 		!Interpolante bilineal
! 		v_g=( f11*(x_2-x0(i))*(y_2-y0(i))+ &
! 		      f21*(x0(i)-x_1)*(y_2-y0(i))+ &
! 		      f12*(x_2-x0(i))*(y0(i)-y_1) + &
! 		      f11*(x0(i)-x_1)*(y0(i)-y_1) )/((x_2-x_1)*(y_2-y_1)) 
! 		     
! 		!------------------z-------------------------------
! 		!Obtener la f de las coordenadas de los vertices
! 		f11=z_global(m1(1)-1,m1(2)-1)
! 		f12=z_global(m1(1)-1,m1(2))
! 		f21=z_global(m1(1)	 ,m1(2)-1)
! 		f22=z_global(m1(1)	 ,m1(2)	)
! 		!Interpolante bilineal
! 		z_g=( f11*(x_2-x0(i))*(y_2-y0(i))+ &
! 		      f21*(x0(i)-x_1)*(y_2-y0(i))+ &
! 		      f12*(x_2-x0(i))*(y0(i)-y_1) + &
! 		      f11*(x0(i)-x_1)*(y0(i)-y_1) )/((x_2-x_1)*(y_2-y_1)) 
! 		
		!5.- Guardar en archivo
		inquire(FILE=fgauges, EXIST=lexist)
		if (.NOT. lexist) then
		  open(20,file=fgauges,status='new',action='write')
! 		  pause
		else    
		  if (it==0.0D0 .and. i==1) then
		   print*,it,i
		   pause
		    open(20,file=fgauges,status='replace',action='write')
		  else
		    open(20,file=fgauges,status='old',action='write', position='append')
		  end if
		end if
		!se podria guardar id,x,y, en otro archivo aparte
		
		write(20,21) id0(i),treal,x0(i),y0(i),z_g,h_g+z_g,u_g,v_g 
		21 format(I	, TR1, F15.5, TR1,	F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5 / )
		close(20)
		
end do
END SUBROUTINE outputgauges


subroutine interpxy(m1,x,y,z,x0,y0,z0)
use global_variables

real (kind=8), dimension(2)::m1 ! par (i,j) del vertice superior derecho de la celda que contiene a x0,y0
real (kind=8)::x0,y0,z0,f11,f21,f12,f22,x_1,x_2,y_1,y_2
real (kind=8),dimension(Nbx,Nby)::x,y,z

!2.-Obtener las coordenadas de los vertices 

x_1=x(m1(1)-1,m1(2))
x_2=x(m1(1),m1(2))
y_1=y(m1(1),m1(2)-1)
y_2=y(m1(1),m1(2))
! pause
!Interpolar!

f11=z(m1(1)-1,   m1(2)-1)
f12=z(m1(1)-1,	  m1(2))
f21=z(m1(1)	 ,m1(2)-1)
f22=z(m1(1)	 ,m1(2)	)
!Interpolante bilineal

z0= ( f11*(x_2-x0)*(y_2-y0)+ &
      f21*(x0-x_1)*(y_2-y0)+ &
      f12*(x_2-x0)*(y0-y_1) + &
      f22*(x0-x_1)*(y0-y_1) )/((x_2-x_1)*(y_2-y_1)) 

end subroutine interpxy
