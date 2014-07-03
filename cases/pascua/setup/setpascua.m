%%
%----------------------------------------------------
%---Parámetros de Discretización y adimensionalizacion
%-----------------------------------------------------
nxi=130;
neta=30;
tinit=0.0;
tfinal=20.0;
cfl=0.95;
caso=999;
batiopt=1;
initqopt=1;
batifilex='data/gridx.dat';
batifiley='data/gridy.dat';
batifilez='data/gridz.dat';
initfileh='data/inith.dat';
initfileu='data/initu.dat';
initfilev='data/initv.dat';
input{1}=num2str(caso);
input{length(input)+1}=[num2str(tinit) 'D0'];%tinit
input{length(input)+1}=[num2str(tfinal) 'D0'];%tfinal
input{length(input)+1}=[num2str(cfl) 'D0' ];%cfl
input{length(input)+1}=num2str(nxi);%nxi
input{length(input)+1}=num2str(neta);%neta
input{length(input)+1}=num2str(batiopt);%batiopt
input{length(input)+1}=batifilex;;%batifile name x
input{length(input)+1}=batifiley;%batifile name y
input{length(input)+1}=batifilez;%batifile name Z
input{length(input)+1}=num2str(initqopt);%batifile name Z
input{length(input)+1}=initfileh;%inith
input{length(input)+1}=initfileu;%initu
input{length(input)+1}=initfilev;%initv
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
dit = 20;
%  dtout=0.1
input{length(input)+1}=num2str(dit);%dit
%  input{length(input)+1}=num2str(dtout)
input{length(input)+1}='1E-10';%kappa, para los ceros numericos
input{length(input)+1}='1';%rk4 
input{length(input)+1}='1';%minmod 
printf('imprimir cada \t dit = %i \t interaciones\n',dit);

%----------------------------------------------------
%--------------------Fricción------------------------
%----------------------------------------------------
input{length(input)+1}='0';%fopt friccion
input{length(input)+1}='1';%outopt 1 = matlab 
outdir='caso1/'
%  outdir=['n',num2str(nxi),'/']
mkdir(outdir)
input{length(input)+1}=outdir;%outdir 
fout=fopen('data/input.dat','w');
for i=1:length(input)
  fprintf(fout,'%s\n', input{i});
end
fclose(fout);
%----------------------------------------------------
%------------------------GEOMETRIA------------------------
%_---------------------------------------------------
%------------------rompimiento de presa-------------

s=load('data/bathy.dat');
x=reshape(s(:,1),neta,nxi);
y=reshape(s(:,2),neta,nxi);
z=reshape(s(:,3),neta,nxi);
h=zeros(size(x));
u=zeros(size(x));
v=zeros(size(x));

h=0.56-z;
h(h<0)=0;

X=X';
Y=Y';
z=z';
h=h';
u=u';
v=v';
save('data/gridx.dat','-ascii','x');
save('data/gridy.dat','-ascii','y');
save('data/gridz.dat','-ascii','z');
save('data/inith.dat','-ascii','h')
save('data/initu.dat','-ascii','u')
save('data/initv.dat','-ascii','v')

save('data/initq.dat','-ascii','huv');
%----------------------------------------------------
%---------------------SERIES DE TIEMPO---------------
%_---------------------------------------------------
b=[1, 5.9, 0.25;
   3, 7.6, 0.25;
   10, 9.644, 0.25;
   15, 10., 0.25;
   22, 10.462, 0.25;
   28, 10.732, 0.25;
   37, 11.005, 0.25;
   38, 11.024, 0.25;
   39, 11.045, 0.25;
   40, 11.12, 0.25;
   46, 11.57, 0.25;]
b(:,2)=b(:,2)+3.;   

fid=fopen('data/gauges.dat','w');
fprintf(fid,'%i\n',size(b,1))
for i=1:size(b,1)
  fprintf(fid,'%i %5.5f %5.5f \n',b(i,1),b(i,2),b(i,3));
end
fclose(fid);
hold on;
figure()
scatter(b(:,2),b(:,3),40)
