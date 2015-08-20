caso=999
tinit=0.0
tfinal=30.0
cfl=0.9
nxi=100
neta=100
batiopt=1
initqopt=1
batifiles=['data/gridX.dat','data/gridY.dat','data/gridZ.dat']
initqfiles=['data/inith.dat', 'data/initu.dat', 'data/initv.dat']
dxi=1.
deta=1.
L=1.
H=1.
U=1.
bcxi0=1
bcxiN=1
bceta0=1
bcetaN=3
dit=10
kappa=1e-5
rktype=1
limtype=1
fricopt=0
outopt=1
outdir='N'+str(nxi)+'P'

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

#-----------geometry and initq-------------
import numpy as np
#import matplotlib.pyplot as plt
x,y=np.meshgrid(np.linspace(-100,100,neta),np.linspace(-100,100,nxi))
z=np.zeros(x.shape)
z[np.nonzero((abs(x)<2.5)* (y <-5))] =20
z[np.nonzero((abs(x)<2.5)* (y >70))] =20
#plt.pcolormesh(x,y,z)

h=5.*np.ones(x.shape)
u=np.zeros(x.shape)
v=np.zeros(x.shape)
h[np.nonzero(x<-2.5)]=10
h[np.nonzero(z>10)]=0

np.savetxt('data/gridX.dat',x)
np.savetxt('data/gridY.dat',y)
np.savetxt('data/gridZ.dat',z)
np.savetxt('data/inith.dat',h)
np.savetxt('data/initu.dat',u)
np.savetxt('data/initv.dat',v)


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
