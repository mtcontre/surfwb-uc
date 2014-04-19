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
nx=908;
ny=656;
tfinal=60*60*3;
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
input{length(input)+1}=num2str(0);%condicion de borde xi_1, 0=custom (soloxi0), 1 = cerrado, 2 = periodic, 3=abierto
%    input{length(input)+1}=num2str(9);%GA 9
%    input{length(input)+1}=num2str(100);%Nsenal
input{length(input)+1}=num2str(3);%condicion de borde xi=nx
input{length(input)+1}=num2str(0);%condicion de borde eta=1
input{length(input)+1}=num2str(0);%condicion de borde eta=ny
printf('borde xi = 1 \t %s \n',input{length(input)-3});
printf('borde xi = nx \t %s\n',input{length(input)-2});
printf('borde eta = 1 \t %s\n',input{length(input)-1});
printf('borde eta = ny \t %s\n',input{length(input)});

%----------------------------------------------------
%----------------------Otros parámetros--------------
%----------------------------------------------------
dit = 5
input{length(input)+1}=num2str(dit);%dit
input{length(input)+1}='1E-6';%kappa, para los ceros numericos
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


