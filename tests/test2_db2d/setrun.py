import numpy as np
import matplotlib.pyplot as plt
import os

os.system('mkdir data')
os.system('mkdir results')
os.system('mkdir results/timeseries')
#=====================================
#computational domain
#=====================================

x = np.linspace(-100,100,50)
y = np.linspace(-100,100,50)
x,y = np.meshgrid(x,y)
z = np.where((y<=-30)*(x<=5)*(x>=-5),15.,0.)
z = np.where((y>=70)*(x<=5)*(x>=-5),15.,z)


h = np.where(x<=0.,10.,5.)
h = np.where(z>0.,0.,h)
u = np.zeros_like(x)
v = np.zeros_like(x)


plt.figure()
plt.pcolormesh(x ,y ,z )
plt.colorbar()
plt.savefig('data/bati.png')

plt.figure()
plt.pcolormesh(x ,y ,h )
plt.colorbar()
plt.savefig('data/inith.png')

plt.figure()
plt.pcolormesh(x ,y ,h +z)
plt.colorbar()
plt.savefig('data/initeta.png')

plt.figure()
plt.pcolormesh(x ,y ,u )
plt.savefig('data/initu.png')
plt.colorbar()

plt.figure()
plt.pcolormesh(x ,y ,v )
plt.savefig('data/initv.png')
plt.colorbar()


np.savetxt('data/gridx.dat',x )
np.savetxt('data/gridy.dat',y )
np.savetxt('data/gridz.dat',z )
np.savetxt('data/inith.dat',h )
np.savetxt('data/initu.dat',u )
np.savetxt('data/initv.dat',v )

#=====================================
#wave gauges
#=====================================
b=[[1,-25., 0.],
   [2,0., 0.],
   [3,25., 0.]]; 
f=open('data/gauges.dat','w')
f.write('%i\n'%len(b))
for i in range(len(b)):
  f.write('%i %.18e %.18e\n'%(b[i][0],b[i][1],b[i][2]))

#=====================================
#input.dat parameters
#=====================================
caso=999
tinit=0.0
tfinal=20.
cfl=0.95
nxi=x.shape[0]
neta=x.shape[1]
batiopt=1
batifiles=['data/gridx.dat','data/gridy.dat','data/gridz.dat']
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
    Nsenal1=txi0.shape[0]
bcxiN=1
bceta0=1
bcetaN=3
dit=-1
dtout = 1.
kappa=1e-3
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
if dit==-1:
  f.write('\n%.3f'%dtout)
f.write('\n%3.20e'%kappa)
f.write('\n%i'%rktype)
f.write('\n%i'%limtype)
f.write('\n%i'%fricopt)
f.write('\n%i'%outopt)
f.write('\n%s\n'%outdir)
f.close()