import numpy as np
import matplotlib.pyplot as plt

def dxvariable(x0,x1,dx0,dxn):
  dxm = 0.5*(dxn+dx0)
  n = int((x1-x0)/dxm)
  i = np.arange(n,dtype=float)
  x = x0 + i*((dxn-dx0)*0.5*(i-1)/(n-1)+dx0)
  return x


rutabase = '../data/bcbati/'
batifiles=[rutabase+'/xmesh.dat',rutabase+'/ymesh.dat',rutabase+'/zmesh.dat']
bcfile = '../data/Seaside_data_final/Seaside_data_V2/01Wavegage.txt'


##bathymetry
xb = np.loadtxt(batifiles[0],unpack=True)
yb = np.loadtxt(batifiles[1],unpack=True)
zb = np.loadtxt(batifiles[2],unpack=True)

dx=0.05
x = dxvariable(xb.min(),30.,8*dx,dx)
x = np.hstack([x,np.arange(30,xb.max(),dx)])

y = dxvariable(yb.min(),-6.,8*dx,dx)
y = np.hstack([y,np.arange(-6.,yb.max(),dx)])

from scipy.interpolate import griddata
x,y = np.meshgrid(x,y,indexing='ij')
z = griddata( np.array([xb.ravel(),yb.ravel()]).T, zb.ravel(), (x,y), method='linear')

##add wall
#dx=0.05
x_extra = np.vstack([x[-1,:]+dx,x[-1,:]+2*dx, x[-1,:]+3*dx, x[-1,:]+4*dx])
y_extra = np.vstack([y[-1,:]   ,y[-1,:]     , y[-1,:]     , y[-1,:]])
z_extra = np.ones(x_extra.shape)*z.max()

x = np.vstack([x,x_extra])
y = np.vstack([y,y_extra])
z = np.vstack([z,z_extra])


batifiles = ['data/gridX.dat','data/gridY.dat','data/gridZ.dat']
np.savetxt(batifiles[0],x)
np.savetxt(batifiles[1],y)
np.savetxt(batifiles[2],z)

#initial condition

h = np.where((z>0.97),0.,0.97-z)
h = np.where(x>35.,0.,h)
u = np.zeros(x.shape)
v = np.zeros(y.shape)

initqfiles=['data/inith.dat', 'data/initu.dat', 'data/initv.dat']
np.savetxt(initqfiles[0],h)
np.savetxt(initqfiles[1],u)
np.savetxt(initqfiles[2],v)

#boundary condition
#50 hz data, dt = 0.02
bc = np.loadtxt(bcfile,skiprows=1)

twg1 = bc[:,0]
etawg1 = np.zeros((x.shape[1],bc.shape[0]))
for it in range(bc.shape[0]):
  etawg1[:,it] = bc[it,3]
np.savetxt('data/etaxi0.dat',etawg1)
np.savetxt('data/timexi0.dat',twg1)
plt.figure()
plt.pcolormesh(etawg1)
plt.savefig('bc.png')
plt.close()

###----------output timeseries-------
b = np.loadtxt('../data/bcbati/outgauges.dat',skiprows=1)
f=open('data/gauges.dat','w')
f.write('%i\n'%len(b))
for i in range(len(b)):
  f.write('%i %.18e %.18e\n'%(int(b[i][0]),b[i][2],b[i][2]))
f.close()



#-----need to check everything
plt.figure(figsize=(8.,4.5))
plt.contour(x,y,z,60,colors='k')
plt.contour(x,y,h,[0.],colors='b')
plt.pcolormesh(x,y,np.ma.masked_where(h<=1e-5,z+h),vmin=0.97-0.3,vmax=0.97+0.3)
plt.plot([30.,30.],[y.min(),y.max()],linewidth=2.,color='k')
plt.colorbar()
plt.scatter(b[:,1],b[:,2],facecolor='r',edgecolor='none')
plt.axis('equal')
#plt.show()
plt.savefig('q0.png')
plt.close()

from mayavi import mlab
mlab.figure()
mlab.mesh(x,y,z)
eta = np.where(h<=1e-5,np.nan,z+h)
mlab.mesh(x,y,eta,opacity=0.4)
#mlab.show()
mlab.savefig('q03d.png',magnification=3)
mlab.close()
plt.figure(figsize=(8.,4.5))
dn = 100
for n in range(0,x.shape[1],dn):
  plt.plot(x[:,n],z[:,n]+h[:,n])
  plt.plot(x[:,n],z[:,n])
plt.tight_layout()
plt.savefig('q0profile.png')


###----------output timeseries-------
b = np.loadtxt('../data/bcbati/outgauges.dat',skiprows=1)
f=open('data/gauges.dat','w')
f.write('%i\n'%len(b))
for i in range(len(b)):
  f.write('%i %.18e %.18e\n'%(int(b[i][0]),b[i][1],b[i][2]))
f.close()

#--------parameters-------
caso=999
tinit=10.5
tfinal=40.
cfl=0.45
nxi=x.shape[0]
neta=x.shape[1]
batiopt=1
initqopt=1
dxi=1.
deta=1.
L=1.
H=1.
U=1.
bcxi0=4
if bcxi0==4:
  GA1=9
  if GA1==9 or GA1==1:
    Nsenal1=twg1.shape[0]
bcxiN=1
bceta0=1
bcetaN=1
dit=-1
if dit==-1:
  dtout=0.01
kappa=1e-5
rktype=1
limtype=1
fricopt=0
outopt=1
outdir='results/'

import os
os.system('mkdir %s'%outdir)
os.system('mkdir %s/timeseries'%outdir)
os.system('mkdir data')
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
  f.write('\n%.18e'%dtout)
f.write('\n%.18e'%kappa)
f.write('\n%i'%rktype)
f.write('\n%i'%limtype)
f.write('\n%i'%fricopt)
f.write('\n%i'%outopt)
f.write('\n%s\n'%outdir)
f.close()
print '---input.dat generated---'


