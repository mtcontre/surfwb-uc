from pylab import *
from colormaps import *
font = {'size'   : 8}
matplotlib.rc('font', **font)
#ion()
p=loadtxt('../results/param.dat')
case=int(p[0])
nx=int(p[1])
ny=int(p[2])
dit=int(p[-1])
t=loadtxt('../results/Time'+str(case)+'.dat')
ntoplot=len(t)-1#min(len(t)-1,10)
rango=range(1,len(t))
x=array([])
y=array([])
figure()
nfr=0
#pt=loadtxt('../data/gauges.dat')
#rango=array([4250,4375,5000,5500,6375,7625,10500,14750])

for i in rango:
  print i
  nfr+=1
  i=int(i)
  
  s=loadtxt('../results/SOL2D.%i.dat.gz' %i)
  if i==rango[0]:
    x=reshape(s[:,0],(ny,nx))
    y=reshape(s[:,1],(ny,nx))
    z=reshape(s[:,2],(ny,nx))
  h=reshape(s[:,3],(ny,nx))
  u=reshape(s[:,4],(ny,nx))
  v=reshape(s[:,5],(ny,nx))
  f=figure()#figsize=(5.5,0.8*3.4)
  pcolormesh(x,y,z,cmap=geo_land)
  pcolormesh(x,y,ma.masked_where(h<=0,z+h),cmap=geo_water,vmin=0.,vmax=10.0)#
  axis('equal')
  colorbar()

  if nfr==10:
    #close('all')
    nfr=0

  f.tight_layout()
  title('T=%.4fs '%(t[i/dit]))#,t[-1]))
  
  savefig('fig'+str(i).rjust(6,'0')+'.png',bbox_inches='tight')
  close()
  

