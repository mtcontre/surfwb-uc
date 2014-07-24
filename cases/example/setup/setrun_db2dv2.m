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
nx=80;
ny=80;
tfinal=80.0;
cfl=0.9;
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
dit = 1
input{length(input)+1}=num2str(dit);%dit
input{length(input)+1}='1E-10';%kappa, para los ceros numericos
input{length(input)+1}='1';%rk4 
input{length(input)+1}='1';%minmod 
printf('imprimir cada \t dit = %i \t interaciones\n',dit);

%----------------------------------------------------
%--------------------Fricción------------------------
%----------------------------------------------------
input{length(input)+1}='0';%fopt friccion
input{length(input)+1}='1';%outopt 1 = matlab 
fout=fopen('../data/input.dat','w');
for i=1:length(input)
  fprintf(fout,'%s\n', input{i});
end
fclose(fout);
%----------------------------------------------------
%------------------------BATI------------------------
%_---------------------------------------------------
%------------------rompimiento de presa-------------
% deben ser 3 columnas x,y,z,  en el archivo data/bathy.dat
%  nx=40;
%  ny=40;
[x,y]=meshgrid(linspace(0,200,nx),linspace(0,200,ny));
z=zeros(size(x));
dx=max(diff(x(1,:)));
z(abs(x-100)<dx & y<95-dx)=20;
z(abs(x-100)<dx & y>170-dx)=20;
%----------------------experimento: rotar los ejes, mantener los puntos
x1=x;
y1=y;
%  x=x*cos(pi/12)-y*sin(pi/12);
%  y=x*sin(pi/12)+y*cos(pi/12);
figure()
pcolor(x,y,z)
xyz=[reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];
save('../data/bathy.dat','-ascii','xyz');
title('Batimetria');
colorbar()
print(gcf,'-dpng','bati.png')
%----------------------------------------------------
%---------------------CONDICIÓN INICIAL--------------
%_---------------------------------------------------
h=5*ones(size(x));
u=zeros(size(y));
v=zeros(size(z));
%  h(x*cos(pi/12)+y*sin(pi/12)<=100 & z==0)=10;
%  h(x*cos(pi/12)+y*sin(pi/12)>100)=5;
h(x<100-dx)=10;
h(abs(x-100)<dx & y<95-dx)=0;
h(abs(x-100)<dx & y>170-dx)=00;
%indices : 22 en adelante queda a la derecha
%indices:  20 para atras quedan a la izquierda		
figure()
pcolor(x,y,h);

title('Condicion inicial para h')
colorbar()
huv=[reshape(h,[],1) reshape(u,[],1) reshape(v,[],1)];
save('../data/initq.dat','-ascii','huv')
print(gcf,'-dpng','h0.png')


%----------------------------------------------------
%---------------------SERIES DE TIEMPO---------------
%_---------------------------------------------------
%Matriz con tres columnas: id  x  y, de los puntos a sacar
b=[1,50,133+45;
   2,50,133;
   3,50,50;
   4,150,133+45;
   5,150,133;
   6,150,50;
   7,50,133+20;
   8,50,133;
   9,50,133-20;
  10,150,133+20;
  11,150,133;
  12,150,133-20];   
   
fid=fopen('../data/gauges.dat','w');
fprintf(fid,'%i\n',size(b,1))
for i=1:size(b,1)
  fprintf(fid,'%i %5.5f %5.5f \n',b(i,1),b(i,2),b(i,3));
end
fclose(fid);
hold on;
scatter(b(:,2),b(:,3),40)