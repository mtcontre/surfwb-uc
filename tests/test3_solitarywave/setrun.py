import numpy as np
import matplotlib.pyplot as plt
import os

os.system('mkdir data')
os.system('mkdir results')
os.system('mkdir results/timeseries')

#=====================================
#incident wave
#=====================================
d0 = 0.2
h0 = 0.07

c = np.sqrt(9.81*(d0+h0))
k = np.sqrt(0.75*h0/d0**3)

t = np.linspace(-3.,6.,200)
x = 0.
eta = h0*(np.cosh(k*(x-c*t)))**(-2)

etaL = np.vstack([t,eta]).T
np.savetxt('data/etaL.dat',etaL)

plt.figure()
plt.plot(etaL[:,0],etaL[:,1])
plt.savefig('data/etaL.png')

#=====================================
#computational domain
#=====================================
n=100.
x=np.linspace(0,10.,n)
dx = np.diff(x)[0]
y=np.linspace(-dx*5,dx*5,10)
x,y=np.meshgrid(x,y,indexing='ij')

z= -d0*np.ones(x.shape)
h = d0*np.ones(x.shape)
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
plt.fill_between(x[:,0],z[:,0],z[:,0]+h[:,0],color='b')
#plt.fill_between(x[:,0],0.*z[:,0],z[:,0],color='k')
plt.tight_layout()
plt.savefig('data/condicion_inicial.png',dpi=300)


#=====================================
#wave gauges
#=====================================
b=[[1,0.,0],\
   [2,3.,0.],\
   [3,6.,0.]]
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
tinit=-3.0
tfinal=20.0
cfl=.45
nxi=x.shape[0]
neta=y.shape[1]
batiopt=1
initqopt=1

dxi=1.
deta=1.
L=1.
H=1.
U=1.
bcxi0=4
if bcxi0==4:
  GA1=1
  if GA1==9 or GA1==1:
    Nsenal1=etaL.shape[0]
bcxiN=1
bceta0=1
bcetaN=1
dit=-1
dtout=0.5
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
