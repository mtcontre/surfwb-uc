import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T = 2.24
a = 1.0
r0 = 0.8
h0 = 0.1
g = 9.81
omega = np.sqrt(8*g*h0/a**2)
#A = (a**2-h0**2)/(a**2+h0**2)
A = (a**4-r0**4)/(a**4+r0**4)

def bati(x,y):
  r = np.sqrt(x**2+y**2)
  return -h0*(1-r**2/a**2)

def analytical(t,x,y):
  r = np.sqrt(x**2+y**2)
  z = -h0*(1-r**2/a**2)
  eta = h0*(np.sqrt(1-A**2)/(1-A*np.cos(omega*t))-1. \
	    -r**2/a**2*((1.-A**2)/(1.-A*np.cos(omega*t))**2-1.))
  u  = 1./(1.-A*np.cos(omega*t))*0.5*omega*x*A*np.sin(omega*t)
  v  = 1./(1.-A*np.cos(omega*t))*0.5*omega*y*A*np.sin(omega*t)
  #eta = np.where(eta<z, z, eta)
  #u = np.where(eta<z, 0., u)
  ##v = np.where(eta<z, 0., v)
  return eta,u,v
def q0(x,y):
  eta,u,v = analytical(0.,x,y)
  return eta,u,v

caso = 999
tinit = 0.0
tfinal = T*20
cfl = 0.7
nxi = 100
neta = 100
batiopt = 1
initqopt = 1
batifiles = ['data/gridX.dat','data/gridY.dat','data/gridZ.dat']
initqfiles = ['data/inith.dat', 'data/initu.dat', 'data/initv.dat']
dxi = 1.
deta = 1.
L = 1.
H = 1.
U = 1.
bcxi0 = 1
bcxiN = 1
bceta0 = 1
bcetaN = 3
dit = 100
kappa = 1e-5
rktype = 1
limtype = 1
fricopt = 0
outopt = 1
outdir = 'N'+str(nxi)+'P'

#---------write to file
f=open('data/input.dat','w')
f.write('%i'%caso)
f.write('\n%.8f %.8f %.8f'%(tinit,tfinal,cfl))
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
f.write('\n%i'%bcxiN)
f.write('\n%i'%bceta0)
f.write('\n%i'%bcetaN)
f.write('\n%i'%dit)
f.write('\n%3.20e'%kappa)
f.write('\n%i'%rktype)
f.write('\n%i'%limtype)
f.write('\n%i'%fricopt)
f.write('\n%i'%outopt)
f.write('\n%s'%outdir)
f.close()


#=====================================
#computational domain
#=====================================
x=np.linspace(-1.3,1.3,nxi)
y=np.linspace(-1.3,1.3,neta)
x,y=np.meshgrid(x,y)
z = bati(x,y)
h,u,v = q0(x,y)
h = np.where(h-z>0,h-z,0.)

np.savetxt(batifiles[0],x)
np.savetxt(batifiles[1],y)
np.savetxt(batifiles[2],z)
np.savetxt(initqfiles[0],h)
np.savetxt(initqfiles[1],u)
np.savetxt(initqfiles[2],v)


#plot stuff
#plt.pcolormesh(x,y,h)
#plt.colorbar()
#plt.show()
#from mpl_toolkits.mplot3d import axes3d
#import matplotlib.pyplot as plt
#from matplotlib import cm
#import numpy as np
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(x,y, h, rstride=10, cstride=10,cmap=cm.coolwarm,shade=True,
#linewidth=0, antialiased=False)
#ax.plot_wireframe(x,y, z, rstride=10, cstride=10,linewidth=0.1,color='green')
#plt.show()

#-----------gauges--------
b=np.array([[1,50,133+20],\
   [2,50,133],\
   [3,50,133-20],\
   [4,150,133+20],\
   [5,150,133],\
   [6,150,133-20]])
b[:,1:3]=b[:,1:3]-100
fid=open('data/gauges.dat','w')
fid.write('%i\n'%b.shape[0])
for i in range(b.shape[0]):
  fid.write('%i %5.5f %5.5f \n'%(b[i,0],b[i,1],b[i,2]))
fid.close()

  

#np.nonzero((abs(x)<2.5) * (y >70) ) = 20
