close all
clear all
caso=999;
input{1}=num2str(caso);
%%
%----------------------------------------------------
%---Parámetros de Discretización y adimensionalizacion
%-----------------------------------------------------
nxi=60;
neta=60;
tinit=0.0;
tfinal=80.0;
cfl=0.9;
batiopt=1;
initqopt=1;
input{length(input)+1}=[num2str(tinit) 'D0 ' num2str(tfinal) 'D0 ' num2str(cfl) 'D0 '];%tfinal
input{length(input)+1}=num2str(nxi);%nxi (filas)
input{length(input)+1}=num2str(neta);%neta (columnas)
input{length(input)+1}=num2str(batiopt); %batiopt
%si 0 lee gridX.dat, gridY.dat y gridZ.dat (binario)
%si 1 lee gridX.dat, gridY.dat y gridZ.dat (ascii)
%si 2 lee un solo archivo con tres columnas: fija *xi* y lee *eta* (hacia la derecha y abajo) (ascii)
%si 3 lee un solo archivo con tres columnas: fija *eta* y lee *xi* (hacia la derecha y abajo) (ascii)
%**si batiopt distinto de 0, lo lee y lo transforma a batiopt=0 
%**guardando archivos gridX.dat, gridY.dat, gridZ.dat, en la carpeta results
%**mejor que sieempre grabe la bati en la carpeta results, pero como  output
if input{end}=='0' || input{end}=='1' 
  input{length(input)+1}='data/gridX.dat';
  input{length(input)+1}='data/gridY.dat';
  input{length(input)+1}='data/gridZ.dat';
elseif input{end}=='2' || input{end}=='3'
  input{length(input)+1}='data/bathy.dat';
end
%lo mismo que batiopt pero para la condicion inicial
%ahora en vez de x,y,z se lee h,u,v
input{length(input)+1}=num2str(initqopt); %initqopt
if input{end}=='0' || input{end}=='1'
  input{length(input)+1}='data/inith.dat';
  input{length(input)+1}='data/initu.dat';
  input{length(input)+1}='data/initv.dat';
elseif input{end}=='2' || input{end}=='3'
  input{length(input)+1}='data/initq.dat';
end 
input{length(input)+1}='1.0D0';%dxi
input{length(input)+1}='1.0D0';%deta
input{length(input)+1}='1.0D0';%L: escala horizontal de longitud
input{length(input)+1}='1.0D0';%H: escala vertical de longitud
input{length(input)+1}='1.0D0';%U: escala de velocidades horizontales

printf('---------init.dat--------\n');
printf('nx = %i \t ny =%i \t nelem = %i \n',neta,nxi,neta*nxi);
printf('t0 = 0.0 \t tfinal = %.2f \t cfl = %.3f\n', tfinal, cfl);

%----------------------------------------------------
%-------Parámetros de condiciones de borde-----------
%----------------------------------------------------
%condiciones de borde implementadas:
%0=custom, 1 = cerrado, 2 = periodic, 3=abierto, 4=GenAbs (?), 5=Inflow-Outflow (?)
%          borde xi=xi(1,:)
%----------------------------------------------------------------------------------------deberia dejar mi condicion de borde como numero 6, no 0..
input{length(input)+1}=num2str(1);%%condicion de borde xi=xi(1,:), 
printf('borde xi = 1 \t %s \n',input{end});
if input{end}=='4'
  %GA=1,2,3,9...??..se que la 9 es la variable y que sí funciona en el borde xi0    
  input{length(input)+1}=num2str(9);%GA
  input{length(input)+1}=num2str(100);%Nsenal 
  %MEJORAS: leer nombre de archivo desde el input.dat------------------------------------------------------
  %...cual nombre de archivo??
end  
%         borde xi=xi(nxi,:)
input{length(input)+1}=num2str(1);%condicion de borde xi=nx
printf('borde xi = nx \t %s\n',input{end});
%         borde eta=eta(:,1)
input{length(input)+1}=num2str(1);%condicion de borde eta=1
printf('borde eta = 1 \t %s\n',input{end});
%         borde eta=eta(:,neta)
input{length(input)+1}=num2str(3);%condicion de borde eta=ny
printf('borde eta = ny \t %s\n',input{end});

%----------------------------------------------------
%----------------------Otros parámetros--------------
%----------------------------------------------------
dit=20
input{length(input)+1}=num2str(dit);%%print every ** iterations
input{length(input)+1}='1E-5';%kappa, para los ceros numericos
input{length(input)+1}='1';%Runge Kutta time stepping method: 1=Rk4, 2=Rk2 
input{length(input)+1}='1';%1=Minmod, 2=Superbee Limiters 
%  printf('imprimir cada \t dit = %i \t in%teraciones\n',dit);

%----------------------------------------------------
%--------------------Fricción------------------------
%----------------------------------------------------
input{length(input)+1}='0';%friction: 0:No, 1:Si
if input{end}=='1'
  input{length(input)+1}='1' %=1: coeficiente unico; =2: matriz de coef: MEJORAR----------------------------
  input{length(input)+1}='1' %tipo de friccion Maning==1, Chezy==2, Sampson==3
  input{length(input)+1}='0.02' %el coeficiente de friccion
end
input{length(input)+1}='1';%formato output 1 = matlab 
fout=fopen('../data/input.dat','w');
for i=1:length(input)
  fprintf(fout,'%s\n', input{i});
end
fclose(fout);
%----------------------------------------------------
%------------------------BATI------------------------
%_---------------------------------------------------
%------------------rompimiento de presa-------------
[x,y]=meshgrid(linspace(-100,100,neta),linspace(-100,100,nxi));
z=zeros(size(x));
dx=max(diff(x(1,:)));
z(abs(x)<dx & y<-5-dx)=20;%-5
z(abs(x)<dx & y>70+dx)=20;%-30
figure()
pcolor(x,y,z)

xyz=[reshape(x',[],1) reshape(y',[],1) reshape(z',[],1)];
save('../data/bathy.dat','-ascii','xyz');%bathi tipo 2

%  xyz=[reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];
%  save('../data/bathy.dat','-ascii','xyz');%bathi tipo 3

save('../data/gridX.dat','-ascii','x');
save('../data/gridY.dat','-ascii','y');
save('../data/gridZ.dat','-ascii','z');


%----------------------------------------------------
%---------------------CONDICIÓN INICIAL--------------
%_---------------------------------------------------
h=5*ones(size(x));
u=zeros(size(y));
v=zeros(size(z));
%  h(x*cos(pi/12)+y*sin(pi/12)<=100 & z==0)=10;
%  h(x*cos(pi/12)+y*sin(pi/12)>100)=5;
h(z>10)=0;
h(x<-dx)=10;
%  h(abs(x)<dx & y<-50-dx)=0;
%  h(abs(x)<dx & y>50+dx)=00;
%  h=0.5*x;
%  u=0.5*y;
%  v=10*x+0.5*y;
figure()
pcolor(x,y,h);
title('Condicion inicial para h')
colorbar()
huv=[reshape(h,[],1) reshape(u,[],1) reshape(v,[],1)];
save('../data/initq.dat','-ascii','huv')%initqopt3

huv=[reshape(h',[],1) reshape(u',[],1) reshape(v',[],1)];
save('../data/initq.dat','-ascii','huv')%initqopt2

save('../data/inith.dat','-ascii','h')%initqopt1
save('../data/initu.dat','-ascii','u')
save('../data/initv.dat','-ascii','v')
%----------------------------------------------------
%---------------------SERIES DE TIEMPO---------------
%_---------------------------------------------------
%Matriz con tres columnas: id  x  y, de los puntos a sacar
b=[1,50,133+20;
   2,50,133;
   3,50,133-20;
   4,150,133+20;
   5,150,133;
   6,150,133-20];
b(:,2:3)=b(:,2:3)-100;
fid=fopen('../data/gauges.dat','w');
fprintf(fid,'%i\n',size(b,1))
for i=1:size(b,1)
  fprintf(fid,'%i %5.5f %5.5f \n',b(i,1),b(i,2),b(i,3));
end
fclose(fid);
hold on;
scatter(b(:,2),b(:,3),40)