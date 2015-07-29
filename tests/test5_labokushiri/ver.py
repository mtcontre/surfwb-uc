import numpy as np
import matplotlib.pyplot as plt
#import colormaps as cmaps
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
os.mkdir('vis')
outdir='results/'
plotdir='vis/'
#read parameters
p=np.loadtxt('%s/param.dat'%outdir)
caso=p[0]
nxi=p[1]
neta=p[2]
nts=p[3]
nit=p[4]
dit=p[5]

#time vector and geometry
t=np.loadtxt('%s/time.dat'%(outdir))
if int(dit)==-1:
  sit=str(1).rjust(8,'0')
  s=np.loadtxt('%s/SOL2D.%s.dat.gz'%(outdir,sit))
else:
  sit=str(0).rjust(8,'0')
  s=np.loadtxt('%s/SOL2D.%s.dat.gz'%(outdir,sit))
x=np.reshape(s[:,0],(nxi,neta),order='F')
y=np.reshape(s[:,1],(nxi,neta),order='F')
z=np.reshape(s[:,2],(nxi,neta),order='F')

#loop
rango=range(0,len(t)-1,1)
#rango = np.nonzero(t>14.)[0]
#rango=[0]
x0=4.7
xf=5.2
y0=1.5
yf=2.2
for i in rango:
  if int(dit)==-1:
    sit=str(i+1).rjust(8,'0')
    s=np.loadtxt('%s/SOL2D.%s.dat.gz'%(outdir,sit))
  else:
    sit=str(i*int(dit)).rjust(8,'0')
    s=np.loadtxt('%s/SOL2D.%s.dat.gz'%(outdir,sit))
    
  h=np.reshape(s[:,3],(nxi,neta),order='F')
  f=plt.figure()#figsize=(5.5,0.8*3.4))
  ax=f.add_subplot(1,1,1)
  p2=ax.pcolormesh(x,y,np.ma.masked_where(h<=1e-5,h+z),vmin=-0.05,vmax=0.05,cmap=plt.cm.coolwarm)
  #ax.set_xticks(np.linspace(y0,yf,5))
  #ax.set_xticklabels(['%.2f'%yi for yi in np.linspace(y0,yf,5)])
  #ax.set_xlim(-yf,-y0)
  #ax.set_ylim(x0,xf)
  ax.set_xlabel('y[m]')
  ax.set_ylabel('x[m]')
  ax.set_aspect(1)
  ax.set_title('T=%ss '%(t[i]))#,t[-1]))
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(p2, cax=cax)
  ax.contour(x,y,z,np.arange(0.,0.2,0.009),colors='k')
  ax.contour(x,y,z,np.arange(-0.04,0.,0.009),colors='k')
  
  f.tight_layout() 
  #show()
  print plotdir+'1dfig'+str(i).rjust(6,'0')+'.png'
  plt.savefig(plotdir+'1dfig'+str(i).rjust(6,'0')+'.png',bbox_inches='tight')
  plt.close()

