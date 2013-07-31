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
nx=115
;
ny=30;
tfinal=20.0;
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
input{length(input)+1}=num2str(0);%condicion de borde xi_1, 0=custom (soloxi0), 1 = cerrado, 2 = periodic, 3=abierto
input{length(input)+1}=num2str(3);%condicion de borde xi=nx
input{length(input)+1}=num2str(1);%condicion de borde eta=1
input{length(input)+1}=num2str(1);%condicion de borde eta=ny
printf('borde xi = 1 \t %s \n',input{length(input)-3});
printf('borde xi = nx \t %s\n',input{length(input)-2});
printf('borde eta = 1 \t %s\n',input{length(input)-1});
printf('borde eta = ny \t %s\n',input{length(input)} );

%----------------------------------------------------
%----------------------Otros parámetros--------------
%----------------------------------------------------
dit = 50
input{length(input)+1}=num2str(dit);%dit
input{length(input)+1}='1E-10';%kappa, para los ceros numericos
input{length(input)+1}='1';%rk4 
input{length(input)+1}='1';%minmod 
printf('imprimir cada \t dit = %i \t interaciones\n',dit);

%----------------------------------------------------
%--------------------Fricción------------------------
%----------------------------------------------------
input{length(input)+1}='0';%fopt friccion
  if input{length(input)} == '1'     
    input{length(input)+1}='1';%fm =1 si se usa un solo coef; = 2 si se usa matriz de coeficientes
    input{length(input)+1}='1';%Cf: Tipo de coeficiente de friccion: 1=manning, 2 =chezy, 3=sampson
    input{length(input)+1}='.02';%coeficiente 
    printf('Se considera friccion\n');
    printf('n = %.6f \n', input{length(input)});
    disp(input{length(input)});
    pause(2)
  end
%tipo de output
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
s=load('SOL2D.0.dat');
x1=(reshape(s(:,1),130,ny));
y1=(reshape(s(:,2),130,ny));
z1=(reshape(s(:,3),130,ny));
x=x1(16:end,:)';
y=y1(16:end,:)';
z=z1(16:end,:)';
%  clear xyz
figure()
plot(x(:,1),y(:,1),'b');hold on;
plot(x(:,end),y(:,end),'r')
plot(x(1,:),y(1,:),'g')
plot(x(end,:),y(end,:),'m')
legend({'xi0','xiN','eta0','etaN'})

xyz=[reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];
save('../data/bathy.dat','-ascii','xyz');
title('Batimetria');
colorbar()
%----------------------------------------------------
%---------------------CONDICIÓN INICIAL--------------
%_---------------------------------------------------
h=zeros(size(x));
u=zeros(size(y));
v=zeros(size(z));
%  
%  
h=0.56-z;
h(x<1.5)=0.85-z(x<1.5);
h(h<=0)=0;
figure()
contourf(x,y,z+h);hold on;
huv=[reshape(h,[],1) reshape(u,[],1) reshape(v,[],1)];
save('../data/initq.dat','-ascii','huv')



%----------------------------------------------------
%---------------------SERIES DE TIEMPO---------------
%_---------------------------------------------------
a=[ 24	2.175	0.72	0.41
    23	2.517	0.79	0.47
    22	3.003	0.81	0.45
    21	3.619	0.78	0.48
    20	4.224	0.83	0.49
    19	4.834	0.95	0.52
    18	5.665	1.04	0.50
    17	6.273	1.03	0.50
    16	6.884	1.10	0.49
    15	7.49	1.50	0.47
    14	8.101	1.92	0.47
    13	8.715	2.04	0.49
    12	9.331	2.05	0.48
    11	9.935	2.21	0.50
    10	10.572	2.42	0.50
    9	11.178	2.68	0.50
    8	11.791	2.66	0.49
    7	12.401	2.44	0.50
    6	13.011	2.38	0.50
    0	13.011	2.38	0.50];
b=[a(end,:);a(1:(end-1),:)];
b=b(:,1:3);
fid=fopen('../data/gauges.dat','w')
for i=1:size(b,1)
  fprintf(fid,'%d %5.5f %5.5f \n',b(i,1),b(i,2),b(i,3));
end
fclose(fid)

scatter(b(:,2),b(:,3),'r')
save('../data/gauges.dat','-ascii','b');