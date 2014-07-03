
close all
clear all
caso=999;
input{1}=num2str(caso);
%%
%----------------------------------------------------
%---Parámetros de Discretización y adimensionalizacion
%-----------------------------------------------------
nx=393+1;
ny=244;
tfinal=80.0;
cfl=0.45;
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
input{length(input)+1}=num2str(4);%condicion de borde xi_1, 0=custom (soloxi0), 1 = cerrado, 2 = periodic, 3=abierto
  input{length(input)+1}=num2str(9);%GA 9
  input{length(input)+1}=num2str(451);%Nsenal
input{length(input)+1}=num2str(1);%condicion de borde xi=nx
input{length(input)+1}=num2str(1);%condicion de borde eta=1
input{length(input)+1}=num2str(1);%condicion de borde eta=ny
printf('borde xi = 1 \t %s \n',input{length(input)-3});
printf('borde xi = nx \t %s\n',input{length(input)-2});
printf('borde eta = 1 \t %s\n',input{length(input)-1});
printf('borde eta = ny \t %s\n',input{length(input)});

%----------------------------------------------------
%----------------------Otros parámetros--------------
%----------------------------------------------------
dit = 25
input{length(input)+1}=num2str(dit);%dit
input{length(input)+1}='1E-5';%kappa, para los ceros numericos
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
s=load('../data/bathy.dat');
x=reshape(s(:,1),ny,nx);
y=reshape(s(:,2),ny,nx);
z=reshape(s(:,3),ny,nx);


%----------------------------------------------------
%---------------------CONDICIÓN INICIAL--------------
%_---------------------------------------------------
q=load('../data/initq.dat');
h=reshape(q(:,1),ny,nx);

	
figure()
pcolor(x,y,h),shading flat,axis equal

figure()
plot(x(:,1),y(:,1),'b');hold on;
plot(x(:,end),y(:,end),'r')
plot(x(1,:),y(1,:),'g')
plot(x(end,:),y(end,:),'m')
legend({'xi0','xiN','eta0','etaN'})
%----------------------------------------------------
%---------------------SERIES DE TIEMPO---------------
%_---------------------------------------------------
%Matriz con tres columnas: id  x  y, de los puntos a sacar
b=[1,4.521, 1.196;
   2,4.521, 1.696;
   3,4.521, 2.196;
   4,0.200, 2.700;
   5,0.200, 1.700;
   6,0.200, 0.500;
   7,3.000, 2.700;
   8,3.000, 1.700;
   9,3.000, 0.500]; 
fid=fopen('../data/gauges.dat','w');
fprintf(fid,'%i\n',size(b,1))
for i=1:size(b,1)
  fprintf(fid,'%i %5.5f %5.5f \n',b(i,1),b(i,2),b(i,3));
end
fclose(fid);
hold on;
%  scatter(b(:,2),b(:,3),40)

