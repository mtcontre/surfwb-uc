from pylab import *
from colormaps import *
#ion()
p=loadtxt('../results/param.dat')
case=int(p[0])
nx=int(p[1])
ny=int(p[2])
dit=int(p[-1])
t=loadtxt('../results/Time'+str(case)+'.dat')
ntoplot=len(t)-1#min(len(t)-1,10)
rango=linspace(0,len(t)*dit,len(t)+1)

x=array([])
y=array([])
figure()
nfr=0
for i in rango:
  print i
  nfr+=1
  i=int(i)
  s=loadtxt('../results/SOL2D.%i.dat.gz' %i)
  if i==0:
    x=reshape(s[:,0],(ny,nx))
    y=reshape(s[:,1],(ny,nx))
    z=reshape(s[:,2],(ny,nx))
  h=reshape(s[:,3],(ny,nx))
  u=reshape(s[:,4],(ny,nx))
  v=reshape(s[:,5],(ny,nx))
  pcolormesh(x,y,z,cmap=geo_land)
  pcolormesh(x,y,ma.masked_where(h<=0,z+h),cmap=geo_water,vmin=-0.03,vmax=0.03)
  colorbar()
  if u.max()>0 or v.max()>0:
    a=1
    #quiver(x,y,ma.masked_where(z>10,u),ma.masked_where(z>10,v))
  if nfr==10:
    #close('all')
    nfr=0
  #axis('tight')
  #xlim(x.min(),x.max())
  #ylim(y.min(),y.max())
  axis('equal')
  title('T=%.4fs de %.4fs'%(t[i/dit],t[-1]))
  #show()
  #if i==0 : colorbar()
  #pause(0.01)
  savefig('fig'+str(i).rjust(6,'0')+'.png')
  close()
  

