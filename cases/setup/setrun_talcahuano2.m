close all
clear all
caso=777;
input{1}=num2str(caso);
%%
%----------------------------------------------------
%---Parámetros de Discretización y adimensionalizacion
%-----------------------------------------------------
nx=600-8;
ny=600-8;
nxi=nx;
neta=ny;
tfinal=130.0*60;
cfl=0.45;
batiopt=1;%matrices ascii
initqopt=1;%matrices ascii
input{length(input)+1}=[num2str(tfinal) 'D0'];%tfinal
input{length(input)+1}=[num2str(cfl) 'D0' ];%cfl
input{length(input)+1}=num2str(nxi);%nxi (filas)
input{length(input)+1}=num2str(neta);%neta (columnas)
input{length(input)+1}=num2str(batiopt); %batiopt
%si 0 lee gridX.dat, gridY.dat y gridZ.dat (binario)
%si 1 lee gridX.dat, gridY.dat y gridZ.dat (ascii)
%si 2 lee un solo archivo con tres columnas: fija *xi* y lee *eta* (hacia la derecha y abajo) (ascii)
%si 3 lee un solo archivo con tres columnas: fija *eta* y lee *xi* (hacia la derecha y abajo) (ascii)
%**si batiopt distinto de 0, lo lee y lo transforma a batiopt=0 
%**guardando archivos gridX.dat, gridY.dat, gridZ.dat, en la carpeta results
%**mejor que sieempre grabe la bati, en binario y en la carpeta results, pero como  output
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
  input{length(input)+1}='data/inith2.dat';
  input{length(input)+1}='data/initu2.dat';
  input{length(input)+1}='data/initv2.dat';
elseif input{end}=='2' || input{end}=='3'
  input{length(input)+1}='data/initq.dat';
end 
input{length(input)+1}='1.0D0';%dxi
input{length(input)+1}='1.0D0';%deta
input{length(input)+1}='1.0D0';%L: escala horizontal de longitud
input{length(input)+1}='1.0D0';%H: escala vertical de longitud
input{length(input)+1}='1.0D0';%U: escala de velocidades horizontales

printf('---------init.dat--------\n');
printf('nx = %i \t ny =%i \t nelem = %i \n',nx,ny,nx*ny);
printf('t0 = 0.0 \t tfinal = %.2f \t cfl = %.3f\n', tfinal, cfl);

%----------------------------------------------------
%-------Parámetros de condiciones de borde-----------
%----------------------------------------------------
%condiciones de borde implementadas:
%0=custom, 1 = cerrado, 2 = periodic, 3=abierto, 4=GenAbs (?), 5=Inflow-Outflow (?)
%          borde xi=xi(1,:)
%----------------------------------------------------------------------------------------deberia dejar mi condicion de borde como numero 6, no 0..
input{length(input)+1}=num2str(0);%%condicion de borde xi=xi(1,:), 
printf('borde xi = 1 \t %s \n',input{end});
if input{end}=='4'
  %GA=1,2,3,9...??..se que la 9 es la variable y que sí funciona en el borde xi0    
  input{length(input)+1}=num2str(9);%GA
  input{length(input)+1}=num2str(100);%Nsenal 
  %MEJORAS: leer nombre de archivo desde el input.dat------------------------------------------------------
  %...cual nombre de archivo??
end  
%         borde xi=xi(nxi,:)
input{length(input)+1}=num2str(0);%condicion de borde xi=nx
printf('borde xi = nx \t %s\n',input{end});
%         borde eta=eta(:,1)
input{length(input)+1}=num2str(0);%condicion de borde eta=1
printf('borde eta = 1 \t %s\n',input{end});
%         borde eta=eta(:,neta)
input{length(input)+1}=num2str(0);%condicion de borde eta=ny
printf('borde eta = ny \t %s\n',input{end});

%----------------------------------------------------
%----------------------Otros parámetros--------------
%----------------------------------------------------
dit=10
input{length(input)+1}=num2str(dit);%%print to file every ** iterations
input{length(input)+1}='1E-5';%kappa, para los ceros numericos
input{length(input)+1}='1';%Runge Kutta time stepping method: 1=Rk4, 2=Rk2 
input{length(input)+1}='1';%1=Minmod, 2=Superbee Limiters 
printf('imprimir cada \t dit = %i \t interaciones\n',dit);

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
%  [x,y]=meshgrid(linspace(0,200,nx),linspace(0,200,ny));
%  z=zeros(size(x));
%  dx=max(diff(x(1,:)));
%  z(abs(x-100)<dx & y<95-dx)=20;
%  z(abs(x-100)<dx & y>170-dx)=20;
%  figure()
%  pcolor(x,y,z)
%  
%  xyz=[reshape(x',[],1) reshape(y',[],1) reshape(z',[],1)];
%  save('../data/bathy.dat','-ascii','xyz');%bathi tipo 2
%  
%  %  xyz=[reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];
%  %  save('../data/bathy.dat','-ascii','xyz');%bathi tipo 3
%  
%  save('../data/gridX.dat','-ascii','x');
%  save('../data/gridY.dat','-ascii','y');
%  save('../data/gridZ.dat','-ascii','z');


%----------------------------------------------------
%---------------------CONDICIÓN INICIAL--------------
%  %_---------------------------------------------------
%  q=load('../results/SOL2D.110960.dat.gz');	
%  h=reshape(q(:,4),nx,ny);
%  u=reshape(q(:,5),nx,ny);
%  v=reshape(q(:,6),nx,ny);
%  %  title('Condiciodn inicial para h')
%  %  colorbar()
%  huv=[reshape(h,[],1) reshape(u,[],1) reshape(v,[],1)];
%  save('../data/initq2.dat','-ascii','huv')%initqopt3
%  
%  huv=[reshape(h',[],1) reshape(u',[],1) reshape(v',[],1)];
%  save('../data/initq2.dat','-ascii','huv')%initqopt2
%  save('../data/inith2.dat','-ascii','h')%initqopt1
%  save('../data/initu2.dat','-ascii','u')
%  save('../data/initv2.dat','-ascii','v')
%----------------------------------------------------
%---------------------SERIES DE TIEMPO---------------
%_---------------------------------------------------
%Matriz con tres columnas: id  x  y, de los puntos a sacar
%  b=[1,665470+(620-1)*6,5932482+(775-1)*6;%mareografo,rll+dr*(nx,ny)
%     2,50,133;
%     3,50,133-20;
%     4,150,133+20;
%     5,150,133;
%     6,150,133-20];
b=[1,668921,5935519;
  2,668909,5935309;
  3,668837,5935047;
  4,668769,5934972;
  5,668580,5934943;
  6,668369,5934934;
  7,668136,5934959;
  8,668127,5935180;
  9,668137,5935405;
 10,668138,5935635;
 11,668155,5935843;
 12,668420,5935901;
 13,668542,5935838;
 14,668389,5935734;
 15,668374,5935531;
 16,668482,5935380;
 17,668614,5935249;
181,668701,5935161;
182,668572,5935084;
 19,668404,5935206;
 20,668308,5935276;
 21,669143,5934940;
 22,669161,5935560;
 23,669136,5936011;
 24,668588,5936036;
 25,668050,5936029];
fid=fopen('../data/gauges.dat','w');
fprintf(fid,'%i\n',size(b,1))
for i=1:size(b,1)
  fprintf(fid,'%i %5.5f %5.5f \n',b(i,1),b(i,2),b(i,3));
end
fclose(fid);
hold on;
scatter(b(:,2),b(:,3),40)
