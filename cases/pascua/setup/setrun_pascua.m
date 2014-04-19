%  tfinal   62.0000000000000     
%   treal+dtreal   62.0102307679407     
%   Iteraciones   5467.00000000000     
%   Time Elapsed =    119.595474000000       seconds.
%  
%  real    2m4.295s
%  user    1m38.670s
%  sys     0m20.989s
%  	

close all
clear all
caso=999;
input{1}=num2str(caso);
%%
%----------------------------------------------------
%---Parámetros de Discretización y adimensionalizacion
%-----------------------------------------------------
nx=130;
ny=30;
tfinal=180;
cfl=0.95;
input{length(input)+1}=[num2str(tfinal) 'D0'];%tfinal	
input{length(input)+1}=[num2str(cfl) 'D0' ];%cfl
input{length(input)+1}=num2str(nx);%nx
input{length(input)+1}=num2str(ny);%ny
input{length(input)+1}='1.0D0';%dxi
input{length(input)+1}='1.0D0';%deta
input{length(input)+1}='1.0D0';%H
input{length(input)+1}='1.0D0';%U	
input{length(input)+1}='1.0D0';%V
printf('---------init.dat--------\n');
printf('nx = %i \t ny =%i \t nelem = %i \n',nx,ny,(nx-1)*(ny-1) );
printf('t0 = 0.0 \t tfinal = %.2f \t cfl = %.3f\n', tfinal, cfl);

%%
%  %----------------------------------------------------
%  %-------Parámetros de condiciones de borde-----------
%  %----------------------------------------------------
input{length(input)+1}=num2str(1);%condicion de borde xi_1, 0=custom (soloxi0), 1 = cerrado, 2 = periodic, 3=abierto
%    input{length(input)+1}=num2str(9);%GA 9
%    input{length(input)+1}=num2str(100);%Nsenal
input{length(input)+1}=num2str(3);%condicion de borde xi=nx
input{length(input)+1}=num2str(1);%condicion de borde eta=1
input{length(input)+1}=num2str(1);%condicion de borde eta=ny
printf('borde xi = 1 \t %s \n',input{length(input)-3});
printf('borde xi = nx \t %s\n',input{length(input)-2});
printf('borde eta = 1 \t %s\n',input{length(input)-1});
printf('borde eta = ny \t %s\n',input{length(input)});

%----------------------------------------------------
%----------------------Otros parámetros--------------
%----------------------------------------------------
dit = 1;
input{length(input)+1}=num2str(dit);%dit
input{length(input)+1}='1E-10';%kappa, para los ceros numericos
input{length(input)+1}='1';%rk4 
input{length(input)+1}='1';%minmod 
printf('imprimir cada \t dit = %i \t interaciones\n',dit);

%----------------------------------------------------
%--------------------Fricción+output------------------------
%----------------------------------------------------
input{length(input)+1}='0';%fopt friccion
input{length(input)+1}='1';%outopt 1 = matlab 
fout=fopen('../data/input.dat','w');
for i=1:length(input)
  fprintf(fout,'%s\n', input{i});
end
fclose(fout);
%----------------------------------------------------
%---------------------SERIES DE TIEMPO---------------
%_---------------------------------------------------
%n total de puntos
%Matriz con tres columnas: id  x  y, de los puntos a sacar
b=[0, 1.00000, 1.00000; 
24, 2.17500, 0.72000 ;
23, 2.51700, 0.79000 ;
22, 3.00300, 0.81000 ;
21, 3.61900, 0.78000 ;
20, 4.22400, 0.83000 ;
19, 4.83400, 0.95000 ;
18, 5.66500, 1.04000 ;
17, 6.27300, 1.03000 ;
16, 6.88400, 1.10000 ;
15, 7.49000, 1.50000 ;
14, 8.10100, 1.92000 ;
13, 8.71500, 2.04000 ;
12, 9.33100, 2.05000 ;
11, 9.93500, 2.21000 ;
10, 10.57200, 2.42000 ;
9, 11.17800, 2.68000 ;
8, 11.79100, 2.66000 ;
7, 12.40100, 2.44000 ;
6, 13.01100, 2.38000 ;
];

s=load('../data/bathy.dat');
x=reshape(s(:,1),ny,nx);
y=reshape(s(:,2),ny,nx);
z=reshape(s(:,3),ny,nx);
%  indices=load('../data/indices.dat');%este lo genere desde fortran
%  for i=1:size(indices,1)
%    ind_base=25+8*(i-1);
%    %parto de 1 celda a la derecha y de ahi contra-reloj
%    ix=indices(i,4)+1;iy=indices(i,5)  ;b=[b;ind_base+1,x(iy,ix),y(iy,ix)];
%    ix=indices(i,4)+1;iy=indices(i,5)+1;b=[b;ind_base+2,x(iy,ix),y(iy,ix)];
%    ix=indices(i,4)  ;iy=indices(i,5)+1;b=[b;ind_base+3,x(iy,ix),y(iy,ix)];
%    ix=indices(i,4)-1;iy=indices(i,5)+1;b=[b;ind_base+4,x(iy,ix),y(iy,ix)];
%    ix=indices(i,4)-1;iy=indices(i,5)  ;b=[b;ind_base+5,x(iy,ix),y(iy,ix)];
%    ix=indices(i,4)-1;iy=indices(i,5)-1;b=[b;ind_base+6,x(iy,ix),y(iy,ix)];
%    ix=indices(i,4)  ;iy=indices(i,5)-1;b=[b;ind_base+7,x(iy,ix),y(iy,ix)];
%    ix=indices(i,4)+1;iy=indices(i,5)-1;b=[b;ind_base+8,x(iy,ix),y(iy,ix)];
%  end
%  fid=fopen('../data/gauges.dat','w');
%  fprintf(fid,'%i\n',size(b,1))
%  for i=1:size(b,1)
%    fprintf(fid,'%i %5.5f %5.5f \n',b(i,1),b(i,2),b(i,3));
%  end
%  fclose(fid);
%  hold on;
%  contourf(x,y,z)
%  scatter(b(:,2),b(:,3))
%  axis equal

