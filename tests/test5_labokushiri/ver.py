import numpy as np
import matplotlib.pyplot as plt
#import colormaps as cmaps
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

os.system('mkdir vis')
outdir='results/'
plotdir='vis/'

#==============================
#Read parameters
#==============================
p=np.loadtxt('%s/param.dat'%outdir)
caso=p[0]
nxi=p[1]
neta=p[2]
nts=p[3]
nit=p[4]
dit=p[5]

#==============================
#time vector and geometry from 
#initial condition
#==============================
t=np.loadtxt('%s/time.dat'%(outdir))
sit=str(0).rjust(8,'0')
s=np.loadtxt('%s/SOL2D.%s.dat.gz'%(outdir,sit))
x=np.reshape(s[:,0],(nxi,neta),order='F')
y=np.reshape(s[:,1],(nxi,neta),order='F')
z=np.reshape(s[:,2],(nxi,neta),order='F')

#==============================
#Loop through iterations
#==============================
rango=range(0,len(t),1)
for i in rango:
  if int(dit)==-1:
    sit=str(i).rjust(8,'0')    
  else:
    sit=str(i*int(dit)).rjust(8,'0')    
  s=np.loadtxt('%s/SOL2D.%s.dat.gz'%(outdir,sit))  
  h=np.reshape(s[:,3],(nxi,neta),order='F')
  f=plt.figure()#figsize=(5.5,0.8*3.4))
  ax=f.add_subplot(1,1,1)
  p2=ax.pcolormesh(x,y,np.ma.masked_where(h<=1e-5,h+z),vmin=-.05,vmax=0.05,cmap=plt.cm.coolwarm)
  ax.set_xlabel('x[m]')
  ax.set_ylabel('y[m]')
  ax.set_aspect(1)
  ax.set_title('T=%ss '%(t[i]))#,t[-1]))
  plt.colorbar(p2)  
  f.tight_layout()
  
  #save images
  print plotdir+'fig'+str(i).rjust(6,'0')+'.png'
  f.savefig(plotdir+'fig'+str(i).rjust(6,'0')+'.png',bbox_inches='tight')
  plt.close()

