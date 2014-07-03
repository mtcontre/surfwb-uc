%   treal+dtreal   62.0102307679407     
%   Iteraciones   5467.00000000000     
%   Time Elapsed =    119.595474000000       seconds.
%  
%  real    2m4.295s
%  user    1m38.670s
%  sys     0m20.989s
%  	


%%
%----------------------------------------------------
%---Parámetros de Discretización y adimensionalizacion
%-----------------------------------------------------
nxi=100;
neta=100;
tinit=0.0;
tfinal=10.0;
cfl=0.9;
caso=999;
input{1}=num2str(caso);
input{length(input)+1}=[num2str(tinit) 'D0'];%tinit
input{length(input)+1}=[num2str(tfinal) 'D0'];%tfinal
input{length(input)+1}=[num2str(cfl) 'D0' ];%cfl
input{length(input)+1}=num2str(nxi);%nxi
input{length(input)+1}=num2str(neta);%neta
input{length(input)+1}='1.0D0';%dxi
input{length(input)+1}='1.0D0';%deta
input{length(input)+1}='1.0D0';%L
input{length(input)+1}='1.0D0';%H	
input{length(input)+1}='1.0D0';%U
printf('---------init.dat--------\n');
printf('nxi = %i \t neta =%i \t nelem = %i \n',nxi,neta,(nxi-1)*(neta-1) );
printf('t0 = 0.0 \t tfinal = %.2f \t cfl = %.3f\n', tfinal, cfl);

%%
%  %----------------------------------------------------
%  %-------Parámetros de condiciones de borde-----------
%  %----------------------------------------------------
input{length(input)+1}=num2str(1);%condicion de borde xi_1, 0=custom (soloxi0), 1 = cerrado, 2 = periodic, 3=abierto
%    input{length(input)+1}=num2str(9);%GA 9
%    input{length(input)+1}=num2str(100);%Nsenal
input{length(input)+1}=num2str(1);%condicion de borde xi=nx
input{length(input)+1}=num2str(1);%condicion de borde eta=1
input{length(input)+1}=num2str(3);%condicion de borde eta=ny
printf('borde xi = 1 \t %s \n',input{length(input)-3});
printf('borde xi = nx \t %s\n',input{length(input)-2});
printf('borde eta = 1 \t %s\n',input{length(input)-1});
printf('borde eta = ny \t %s\n',input{length(input)});

%----------------------------------------------------
%----------------------Otros parámetros--------------
%----------------------------------------------------
dit = -1
dtout=0.05
input{length(input)+1}=num2str(dit);%dit
input{length(input)+1}=num2str(dtout)
input{length(input)+1}='1E-10';%kappa, para los ceros numericos
input{length(input)+1}='1';%rk4 
input{length(input)+1}='1';%minmod 
printf('imprimir cada \t dit = %i \t interaciones\n',dit);

%----------------------------------------------------
%--------------------Fricción------------------------
%----------------------------------------------------
input{length(input)+1}='0';%fopt friccion
input{length(input)+1}='1';%outopt 1 = matlab 
outdir=['n',num2str(nxi),'/']
mkdir(['../' outdir])
input{length(input)+1}=outdir;%outdir 
fout=fopen('../data/input.dat','w');
for i=1:length(input)
  fprintf(fout,'%s\n', input{i});
end
fclose(fout);
%----------------------------------------------------
%------------------------GEOMETRIA------------------------
%_---------------------------------------------------
%------------------rompimiento de presa-------------

xi=linspace(-100,100,nxi);
eta=linspace(-100,100,neta);
[x,y]=meshgrid(eta,xi);
z=zeros(size(x));
z(abs(x)<=4 & y >=70)=20;
z(abs(x)<=4 & y <=-5)=20;

h=10*ones(size(x));
u=zeros(size(y));
v=zeros(size(z));
h(x>=-4)=5;
h(z>0)=0.;

figure()
pcolor(x,y,z)
colorbar()
figure()
pcolor(x,y,h);
colorbar()
xyz=[reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];
huv=[reshape(h,[],1) reshape(u,[],1) reshape(v,[],1)];
save('../data/bathy.dat','-ascii','xyz');
save('../data/initq.dat','-ascii','huv')
