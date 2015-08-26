import numpy as np
import matplotlib.pyplot as plt
from string import split
import os

outdir = 'results/'
plt.plotdir = 'ver/'

if not os.path.isdir(plt.plotdir):
  print '%s no existe. Creando directorio'%plt.plotdir
  os.system('mkdir %s'%plt.plotdir)
  
# Read general parameters
f=open(outdir+'grids/gridproperties.dat')

line=f.readline();line=split(line)
dit=int(line[1])

line=f.readline();line=split(line)
nproc=int(line[1])

line=f.readline();line=split(line)
dims=[0,0]
dims[0]=int(line[1])
dims[1]=int(line[2])

line=f.readline();line=split(line)
ngrids=int(line[1])

line=f.readline();line=split(line)
nxi=int(line[1])

line=f.readline();line=split(line)
neta=int(line[1])

batinames = [0,0,0]
line=f.readline();line=split(line);
print line
line = split(line[0],sep='/')
batinames = line[0]
#line=f.readline();line=split(line);line = split(line[0],sep='/')
#batinames[1] = line[1]
#line=f.readline();line=split(line);line = split(line[0],sep='/')
#batinames[2] = line[1]
f.close()


#set the range of frames to plt.plot
t=np.loadtxt(outdir+'time.dat')
rango=range(len(t)-1)


for it in rango:
  #get bathymetry filenames  for each level
  x=np.loadtxt(outdir+''+batinames[0])
  y=np.loadtxt(outdir+''+batinames[1])
  z=np.loadtxt(outdir+''+batinames[2])
  for i in range(dims[0]):
    for j in range(dims[1]):
      #file with properties for this grid
      fname=outdir+'grids/grid%03d_%03d.dat'%(i,j)
      f=open(fname)
      line=f.readline()      
      nbx=int(split(line)[0]);line=f.readline()
      nby=int(split(line)[0]);line=f.readline()
      #print nbx,nby
      si=int(split(line)[0])-1;line=f.readline()#minus one for python
      ei=int(split(line)[0]);line=f.readline()
      sj=int(split(line)[0])-1;line=f.readline()
      ej=int(split(line)[0]);line=f.readline()
      f.close()
      
      fname ='%s/SOL2D.P%03d_%03d_%03d.%08d.dat.gz'%(outdir,nproc,i,j,it)
      #fname=outdir+'P%03dframe%08d.%03d_%03d.dat'%(it,i,j)
      s=np.loadtxt(fname)
      h=np.reshape(s[:,0],(nbx,nby),order='F')
      u=np.reshape(s[:,1],(nbx,nby),order='F')
      v=np.reshape(s[:,2],(nbx,nby),order='F')

      plt.pcolormesh(x[si:(ei+1),sj:(ej+1)],y[si:(ei+1),sj:(ej+1)],h,vmin=-0.08,vmax=0.08)

      if i == dims[0]-1:
	ei-=1
      if j == dims[1] -1:
	ej-=1

      plt.plot(x[si:ei,sj],y[si:ei,sj],color='k')#,linewidth=1.)
      plt.plot(x[si:ei,ej],y[si:ei,ej],color='k')#,linewidth=1.)
      plt.plot(x[si,sj:ej],y[si,sj:ej],color='k')#,linewidth=1.)
      plt.plot(x[ei,sj:ej],y[ei,sj:ej],color='k')#,linewidth=1.)
  plt.axis('equal')
  plt.colorbar()  
  plt.title('P=%i, t=%.3fs'%(nproc,t[it]))
  plt.savefig(plt.plotdir+'/p%03dframexy%08d.png'%(nproc,it))
  plt.close()


#axis('equal')
#
#close()
#if __name__=='__main__':
  #doplt.plots('results/','ver/')
