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
integer, dimension(2) :: m1_temp
write(ncaso,*) caso
ncaso=adjustl(ncaso)
numbercaso=ncaso  
path=outdir
write(ncaso,*) caso
ncaso=(adjustl(ncaso))
fgauges=trim(path)//'/timeseries/gauges.' // trim(ncaso)// '.dat'


do i=1,Nts

  !3.- Generar los datos interpolados (interpolacion bilineal)
  !------------------h-------------------------------
  m1_temp = m1(:,i)
  call interpxy(m1_temp,x_global,y_global,qreal_global(1,:,:),x0(i),y0(i),h_g,0)
  call interpxy(m1_temp,x_global,y_global,qreal_global(2,:,:),x0(i),y0(i),u_g,0)
  call interpxy(m1_temp,x_global,y_global,qreal_global(3,:,:),x0(i),y0(i),v_g,0)
  call interpxy(m1_temp,x_global,y_global,z_global(:,:),x0(i),y0(i),z_g,0)

  !5.- Guardar en archivo
  inquire(FILE=fgauges, EXIST=lexist)
  if (.NOT. lexist) then
    open(20,file=fgauges,status='new',action='write')
  else    
    if (it==0 .and. i==1) then
      open(20,file=fgauges,status='replace',action='write')
    else
      open(20,file=fgauges,status='old',action='write', position='append')
    end if
  end if
  !se podria guardar id,x,y, en otro archivo aparte
  
  write(20,21) id0(i),treal,x0(i),y0(i),z_g,h_g+z_g,u_g,v_g 
  21 format(I3.3,TR1,F15.5, TR1,F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5 / )
  close(20)

end do

END SUBROUTINE outputgauges


subroutine interpxy(m1,x,y,z,x0,y0,z0,order)
use global_variables
integer ::order !que orden de innterpolacion, 0 o 1 
integer, dimension(2)::m1 ! par (i,j) del vertice superior derecho de la celda que contiene a x0,y0
real (kind=8)::x0,y0,z0,f11,f21,f12,f22,x_1,x_2,y_1,y_2
real (kind=8),dimension(Nbx,Nby)::x,y,z

!2.-Obtener las coordenadas de los vertices 

if (order.eq.1) then
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
elseif (order.eq.0) then
  z0=z(m1(1),m1(2))
end if
end subroutine interpxy
