import numpy as np
import matplotlib.pyplot as plt
from string import split
import os

def get_gridprops(outdir):  
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
  line=f.readline(); line=split(line);
  batinames[0] = split(line[0],sep='/')[2]
  batinames[1] = split(line[1],sep='/')[2]
  batinames[2] = split(line[2],sep='/')[2]
  f.close()
  
  return dit,nproc,dims,ngrids,nxi,neta,batinames
def plotframes(outdir,plotdir):

  if not os.path.isdir(plotdir):
    print '%s no existe. Creando directorio'%plotdir
    os.system('mkdir %s'%plotdir)
    
  # Read general parameters
  dit,nproc,dims,ngrids,nxi,neta,batinames = get_gridprops(outdir)


  #set the range of frames to plt.plot
  t=np.loadtxt(outdir+'timeP%03d.dat'%nproc)
  rango=range(len(t))


  for it in rango:
    print 'frame %i/%i of sim. with %i procs and grid %i'%(it,rango[-1],nproc,nxi)
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
	
	fname ='%s/P%03dframe%08d.%03d_%03d.dat'%(outdir,nproc,it,i,j)
	#fname=outdir+'P%03dframe%08d.%03d_%03d.dat'%(it,i,j)
	s=np.loadtxt(fname)
	h=np.reshape(s[:,0],(nbx,nby),order='F')
	u=np.reshape(s[:,1],(nbx,nby),order='F')
	v=np.reshape(s[:,2],(nbx,nby),order='F')

	eta = np.ma.masked_where(h<=1e-3,h+z[si:ei,sj:ej])
	plt.pcolormesh(x[si:(ei+1),sj:(ej+1)],y[si:(ei+1),sj:(ej+1)],
		      eta,vmin=-0.08,vmax=0.05)

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
    plt.savefig(plotdir+'/p%03dframexy%08d.png'%(nproc,it))
    plt.close()


#axis('equal')
#
#close()
if __name__=='__main__':
  import sys
  import os
  arg = sys.argv[1:]
  n = int(arg[0])
  p = int(arg[1])
  print 'to plot with n=%i and p=%i'%(n,p)
  plotframes('results_%i_1_%i/'%(n,p),'ver_%i/'%(n))
