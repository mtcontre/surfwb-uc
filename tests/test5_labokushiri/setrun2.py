import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import urllib
os.system('mkdir data')
os.system('mkdir results')
os.system('mkdir results/timeseries')

#====================================
#get original data from noaa website
#=====================================

rooturl='http://nctr.pmel.noaa.gov/benchmark/Laboratory/Laboratory_MonaiValley/'
batifname0 = 'MonaiValley_Bathymetry.txt'
wavefname0 = 'MonaiValley_InputWave.txt'
if not os.path.isfile('data/'+batifname0):
  urllib.urlretrieve(rooturl+batifname0,'data/'+batifname0)
if not os.path.isfile('data/'+wavefname0):
  urllib.urlretrieve(rooturl+wavefname0,'data/'+wavefname0)
bati = np.loadtxt('data/MonaiValley_Bathymetry.txt',skiprows=1)
wave = np.loadtxt('data/MonaiValley_InputWave.txt',skiprows=1)
#=====================================
#computational domain
#=====================================

ny = np.unique(bati[:,1]).shape[0]
x = np.reshape(bati[:,0],(-1,ny))[::4,::4]
y = np.reshape(bati[:,1],(-1,ny))[::4,::4]
z = -1.*np.reshape(bati[:,2],(-1,ny))[::4,::4]
dx = np.diff(x[:,0])[0]

xwall = x[:,0].max()*np.ones((2,x.shape[1]))
xwall[0,:] += dx
xwall[1,:] *= 2*dx

ywall = np.zeros_like(xwall)
ywall[0,:] = y[-1,:]
ywall[1,:] = y[-1,:]
zwall = np.ones_like(xwall)*0.5

x = np.vstack([x,xwall])
y = np.vstack([y,ywall])
z = np.vstack([z,zwall])


h = np.where(z<=0.,-z,0.)
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
b=[[1,4.521, 1.196],
   [2,4.521, 1.696],
   [3,4.521, 2.196],
   [4,0.000, 2.700],
   [5,0.000, 1.700],
   [6,0.000, 0.500],
   [7,3.000, 2.700],
   [8,3.000, 1.700],
   [9,3.000, 0.500]]; 
f=open('data/gauges.dat','w')
f.write('%i\n'%len(b))
for i in range(len(b)):
  f.write('%i %.18e %.18e\n'%(b[i][0],b[i][1],b[i][2]))

#=====================================
#input.dat parameters
#=====================================
caso=999
tinit=0.0
tfinal=50.
cfl=0.45
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
bcetaN=1
dit=-1
dtout = 0.5
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