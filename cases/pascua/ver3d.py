from pylab import *

#from cmap import *
from mayavi import mlab
#obj.scene.disable_render = True
#mlab.options
#mlab.options.offscreen = True
#ion()
p=loadtxt('../results/param.dat')
case=int(p[0])
nx=int(p[1])
ny=int(p[2])
dit=int(p[-1])
t=loadtxt('../results/Time'+str(case)+'.dat')
ntoplot=len(t)-1#min(len(t)-1,10)
#rango=linspace(0,ntoplot*dit,ntoplot+1)
rango=arange(0,ntoplot,5*dit)

x=array([])
y=array([])
#figure()
nfr=0
for i in rango:
  #i=100
  #print i
  nfr+=1
  i=int(i)
  s=loadtxt('../results/SOL2D.%i.dat' %i)
  print '%i'%i
  x=reshape(s[:,0],(ny,nx))
  y=reshape(s[:,1],(ny,nx))
  z=reshape(s[:,2],(ny,nx))
  h=reshape(s[:,3],(ny,nx))
  u=reshape(s[:,4],(ny,nx))
  v=reshape(s[:,5],(ny,nx))
  water=where(h>1.e-3,z+h,nan)
  if True:
    mlab.figure(size=(800*0.8,700*0.8))
    mesh_water=mlab.mesh(x.transpose(), y.transpose(),water.transpose(),colormap='Blues',vmin=0.0,vmax=10.1)
    mesh_land=mlab.mesh(x.transpose(),y.transpose(),z.transpose(),colormap='gist_earth',vmin=0.,vmax=1.)
    # Retrieve the LUT of the surf object.
    #lut = mesh_water.module_manager.scalar_lut_manager.lut.table.to_array()
    #lut[:, 0] = geo_water[0]
    #lut[:, 1] = geo_water[1]
    #lut[:, 2] = geo_water[2]
    #mesh_water.module_manager.scalar_lut_manager.lut.table = lut
    #mlab.draw()
    #mlab.colorbar()
    
    #lut = mesh_land.module_manager.scalar_lut_manager.lut.table.to_array()
    #lut[:, 0] = geo_land[0]
    #lut[:, 1] = geo_land[1]
    #lut[:, 2] = geo_land[2]
    #lut[:, 3]= np.linspace(255, 0, 256)
    
    #mesh_land.module_manager.scalar_lut_manager.lut.table = lut
    #mesh_water.actor.actor.scale=[1,1,10]
    #mesh_land.actor.actor.scale=[1,1,7]
    #mlab.draw()
    mlab.colorbar()
    mlab.view(44.702487859422739,
    57.422360898293746,
    662.80712520602356,
    array([  72.22287172,  101.25175099,   41.53355985]))
    mlab.title('T=%.4fs de %.4fs'%(t[i/dit],t[-1]),size=0.5)
    mlab.savefig('fig'+str(i).rjust(6,'0')+'.png')
    mlab.close()

   #su.mlab_source.scalars =h.transpose()
  
  #axis('tight')
  #xlim(x.min(),x.max())
  #ylim(y.min(),y.max())

  
  
  #mlab.axes(extent=[x.min(),x.max(),x.min(),x.max(),9,11])

 
 #709.66679368099665,
 #array([ 100.,  100.,   50.]))

  #mlab.view(distance=800)
  #pause(0.01)

  
  #close()
  

