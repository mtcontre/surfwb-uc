#--------parameters-------
caso=999
tinit=0.0
tfinal=50.0
cfl=0.45
nxi=393+1
neta=244
batiopt=1
batifiles=['data/gridX.dat','data/gridY.dat','data/gridZ.dat']
initqopt=1
initqfiles=['data/inith.dat', 'data/initu.dat', 'data/initv.dat']
dxi=1.
deta=1.
L=1.
H=1.
U=1.
bcxi0=4
if bcxi0==4:
  GA1=9
  if GA1==9 or GA1==1:
    Nsenal1=451
bcxiN=1
bceta0=1
bcetaN=1
dit=25
kappa=1e-5
rktype=1
limtype=1
fricopt=0
outopt=1
outdir='results/'

#---------write to file
f=open('../data/input.dat','w')
f.write('%i'%caso)
f.write('\n%.8f'%tinit)
f.write('\n%.8f'%tfinal)
f.write('\n%.8f'%cfl)
f.write('\n%i'%nxi)
f.write('\n%i'%neta)
f.write('\n%i'%batiopt)
for i in range(len(batifiles)):
  f.write('\n%s'%batifiles[i])
f.write('\n%i'%initqopt)
for i in range(len(initqfiles)):
  f.write('\n%s'%initqfiles[i])
f.write('\n%.8f'%dxi)
f.write('\n%.8f'%deta)
f.write('\n%.8f'%L)
f.write('\n%.8f'%H)
f.write('\n%.8f'%U)
f.write('\n%i'%bcxi0)
if bcxi0==4:
  f.write('\n%i'%GA1)
  if GA1==9 or GA1==1:
    f.write('\n%i'%Nsenal1)
f.write('\n%i'%bcxiN)
f.write('\n%i'%bceta0)
f.write('\n%i'%bcetaN)
f.write('\n%i'%dit)
f.write('\n%3.20e'%kappa)
f.write('\n%i'%rktype)
f.write('\n%i'%limtype)
f.write('\n%i'%fricopt)
f.write('\n%i'%outopt)
f.write('\n%s\n'%outdir)
f.close()


#---------bati+q0 files------
s=np.loadtxt('../data/bathy.dat')
x=np.reshape(s[:,0],xi,neta)

