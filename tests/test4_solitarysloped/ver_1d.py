import numpy as np
import matplotlib.pyplot as plt
import os
os.system('mkdir vis')
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
#rango=[0]
for i in rango:
  if int(dit)==-1:
    sit=str(i).rjust(8,'0')
  else:
    sit=str(i*int(dit)).rjust(8,'0')
  s=np.loadtxt('%s/SOL2D.%s.dat.gz'%(outdir,sit))
  h=np.reshape(s[:,3],(nxi,neta),order='F')
  #u=reshape(s[:,4],(ny,nx))
  #v=reshape(s[:,5],(ny,nx))
  f=plt.figure(figsize=(8.,3.))#figsize=(5.5,0.8*3.4))
  #plot(x[1,:],ma.masked_where(h[1,:]<=0,z[1,:]+h[1,:]))
  plt.fill_between(x[:,1],z.min()*np.ones_like(z[:,1]),z[:,1],color='k')
  plt.fill_between(x[:,1],z[:,1],z[:,1]+h[:,1],color='b')
  
  #plt.fill_between(np.array([x.max(),x.max()+0.1]),np.array([0.,0.]),[0.3,0.3],color=3*[.5])
  #plt.fill_between(np.array([x.min(),x.min()-0.1]),np.array([0.,0.]),[0.3,0.3],color=3*[.5])
  #plt.axis('equal')
  #plt.xlim(-0.15,10.15)
  plt.ylim(-.2,0.3)
  plt.title('T=%ss '%(t[i]))#,t[-1]))
  f.tight_layout() 
  #show()
  print plotdir+'1dfig'+str(i).rjust(6,'0')+'.png'
  plt.savefig(plotdir+'1dfig'+str(i).rjust(6,'0')+'.png',bbox_inches='tight')
  plt.close()

