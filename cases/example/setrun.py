#--------parameters-------
caso=999
tinit=0.0
tfinal=10.0
cfl=0.95
nxi=100
neta=100
batiopt=1
batifiles=['data/gridX.dat','data/gridY.dat','data/gridZ.dat']
initqopt=1
initqfiles=['data/inith.dat', 'data/initu.dat', 'data/initv.dat']
dxi=1.
deta=1.
L=1.
H=1.
U=1.
bcxi0=1
if bcxi0==4:
  GA1=9
  if GA1==9 or GA1==1:
    Nsenal1=451
bcxiN=1
bceta0=1
bcetaN=3
dit=10
kappa=1e-5
rktype=1
limtype=1
fricopt=0
outopt=1
outdir='results/'

#---------write to file
f=open('data/input.dat','w')
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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
x=np.linspace(-100,100,nxi)
y=np.linspace(-100,100,neta)
x,y=np.meshgrid(x,y)
z = np.zeros(x.shape)
z = np.where( (np.abs(x)<2.5)*(y<-70.), 15., z)
z = np.where( (np.abs(x)<2.5)*(y>5.), 15.,z)
h = 5.*np.ones(x.shape)
u = np.zeros(x.shape)
v = np.zeros(x.shape)

h = np.where( x<-2.5,10.,h)

h = np.where( z>0.,0.,h)

plt.pcolormesh(x,y,z)
plt.colorbar()
plt.savefig('data/bati.png')
plt.close()

np.savetxt(batifiles[0],x)
np.savetxt(batifiles[1],y)
np.savetxt(batifiles[2],z)
np.savetxt(initqfiles[0],h)
np.savetxt(initqfiles[1],u)
np.savetxt(initqfiles[2],v)


##----------series de tiempo-------
b=[[1,0.,50.],\
   [2,50.,50.]];
f=open('data/gauges.dat','w')
f.write('%i\n'%len(b))
for i in range(len(b)):
  f.write('%i %.18e %.18e\n'%(b[i][0],b[i][2],b[i][2]))
f.close()