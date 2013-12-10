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
rango=linspace(0,len(t)*dit,len(t)+1)
rango=range(0,len(t)*dit,dit*5)
x=array([])
y=array([])
figure()
nfr=0
pt=loadtxt('../data/gauges.dat')
rango=array([4250,4375,5000,5500,6375,7625,10500,14750])

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
  f=figure(figsize=(5.5,0.8*3.4))
  pcolormesh(x,y,z,cmap=geo_land)
  pcolormesh(x,y,ma.masked_where(h<=0,z+h),cmap=geo_water,vmin=-0.03,vmax=0.03)
  colorbar()
  for k in range(pt.shape[0]):
    xn=pt[k][1]
    yn=pt[k][2]
    plot(xn,yn,'ko',markersize=5)
    text(xn-0.1,yn+0.1,'%i' %pt[k][0],fontsize=10)
    #text(xn+0.1,yn+0.,'(%.2f,%.2f)' %(pt[i][1],pt[i][2]),fontsize=6)
  #clabel(cs,inline=2,fontsize=8)
  if nfr==10:
    #close('all')
    nfr=0
  #axis('tight')
  #xlim(x.min(),x.max())s
  #ylim(y.min(),y.max())
  f.tight_layout() 
  axis('tight')
  f.tight_layout()
  title('T=%.4fs '%(t[i/dit]))#,t[-1]))
  #show()
  #if i==0 : colorbar()
  #pause(0.01)
  
  savefig('00fig'+str(i).rjust(6,'0')+'.png',bbox_inches='tight')
  close()
  

