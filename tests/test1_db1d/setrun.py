import numpy as np
import matplotlib.pyplot as plt
import os

os.system('mkdir data')
os.system('mkdir results')
os.system('mkdir results/timeseries')

#=====================================
#computational domain
#=====================================
dx=0.45/20
n=int(5.6/dx)

x=np.linspace(0,5.6,n)
y=np.linspace(-dx,dx,3)
x,y=np.meshgrid(x,y)

z=np.zeros(x.shape)
z = np.where( (x>=2.39+1.61)*(x<=2.39+1.61+0.45),0.065/0.45*(x-2.39-1.61),z)
z = np.where( (x>=2.39+1.61+0.45)*(x<=2.39+1.61+0.45*2),-0.065/0.45*(x-2.39-1.61-0.45)+0.065,z)

h = np.zeros(x.shape)
h = np.where(x<=2.39,0.111,h)
h = np.where((x>=2.39+1.61+0.45),0.02-z,h)
h = np.where(h<=0,0.,h)

u = np.zeros(x.shape)
v = np.zeros(x.shape)

batifiles=['data/gridX.dat','data/gridY.dat','data/gridZ.dat']
initqfiles=['data/inith.dat', 'data/initu.dat', 'data/initv.dat']

np.savetxt(batifiles[0],x)
np.savetxt(batifiles[1],y)
np.savetxt(batifiles[2],z)
np.savetxt(initqfiles[0],h)
np.savetxt(initqfiles[1],u)
np.savetxt(initqfiles[2],v)

plt.figure(figsize=(8.,3.))
plt.fill_between(x[0,:],z[0,:],z[0,:]+h[0,:],color='b')
plt.fill_between(x[0,:],0.*z[0,:],z[0,:],color='k')
plt.tight_layout()
plt.savefig('data/condicion_inicial.png',dpi=300)


#=====================================
#wave gauges
#=====================================
b=[[1,5.575,0],\
   [2,4.925,0.],\
   [3,3.935,0.]]
f=open('data/gauges.dat','w')
f.write('%i\n'%len(b))
for i in range(len(b)):
  s='%i %.18e %.18e\n'%(b[i][0],b[i][1],b[i][2])
  f.write(s)
f.close()

#=====================================
#input.dat parameters
#=====================================
caso=999
tinit=0.0
tfinal=50.0
cfl=1.
nxi=x.shape[0]
neta=y.shape[1]
batiopt=1
initqopt=1

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
bcetaN=1
dit=-1
dtout=1.
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
if dit==-1:
  f.write('\n%.3f'%dtout)
f.write('\n%3.20e'%kappa)
f.write('\n%i'%rktype)
f.write('\n%i'%limtype)
f.write('\n%i'%fricopt)
f.write('\n%i'%outopt)
f.write('\n%s\n'%outdir)
f.close()
