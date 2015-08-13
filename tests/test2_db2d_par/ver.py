from pylab import *
from string import split

#def doplots(outdir,plotdir):
if True:
  (outdir,plotdir) = ('results/','VN100P1/')
  import os
  if not os.path.isdir(plotdir):
    print '%s no existe. Creando directorio'%plotdir
    os.system('mkdir %s'%plotdir)
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
  nxi=int(line[1])

  line=f.readline();line=split(line)
  neta=int(line[1])

  #set the range of frames to plot
  l=1
  t=loadtxt(outdir+'time.dat')
  frames=range(0,len(t)*dit,dit)
  rango=arange(0,len(frames),l)

  #batinames=[]
  #for lgrid in range(ngrids):
    #line=f.readline();line=split(line)
    #batinames.extend([line])
    #batinames[lgrid]=[split(batinames[lgrid][i],sep='/')[1] for i in range(3)]
  #f.close()

  for i in rango:
    if i<=len(frames)-1:
      it=frames[i]
    else:
      break

    #get bathymetry filenames  for each level
    #for l in range(ngrids):
      #x=loadtxt(outdir+''+batinames[l][0])
      #y=loadtxt(outdir+''+batinames[l][1])
      #z=loadtxt(outdir+''+batinames[l][2])
    for i in range(dims[0]):
      for j in range(dims[1]):
	#file with properties for this grid
	fname=outdir+'grids/grid%03d_%03d.dat'%(i,j)
	f=open(fname)
	line=f.readline()
	nbx=int(split(line)[0]);line=f.readline()
	nby=int(split(line)[0]);line=f.readline()
	#si=int(split(line)[0])-1;line=f.readline()#minus one for python
	#ei=int(split(line)[0]);line=f.readline()
	#sj=int(split(line)[0])-1;line=f.readline()
	#ej=int(split(line)[0]);line=f.readline()
	f.close()
	
	fname ='%s/SOL2D.%03d_%03d.%08d.dat.gz'%(outdir,i,j,it)
	#fname=outdir+'P%03dframe%08d.%03d_%03d.dat'%(it,i,j)
	s=loadtxt(fname)
	x=reshape(s[:,0],(nbx,nby))
	y=reshape(s[:,1],(nbx,nby))
	z=reshape(s[:,2],(nbx,nby))	  
	h=reshape(s[:,3],(nbx,nby))
	u=reshape(s[:,4],(nbx,nby))
	v=reshape(s[:,5],(nbx,nby))

	eta = ma.masked_where(h<=1e-2,z+h)
	pcolormesh(x,y,eta,vmin=0.,vmax=10.)

	#ei=ei+1

	#plot(x[si:ei,sj],y[si:ei,sj],color='k',linewidth=1.)
	#plot(x[si:ei,ej],y[si:ei,ej],color='k',linewidth=1.)
	#plot(x[si,sj:ej],y[si,sj:ej],color='k',linewidth=1.)
	#plot(x[ei,sj:ej],y[ei,sj:ej],color='k',linewidth=1.)
	#axis('equal')
    plt.colorbar()
    title('t=%.3fs'%(t[it/dit]))
    axis('equal')
    savefig(plotdir+'/framexy%08d.png'%it)
    close()
if __name__=='__main__':
  doplots('results/','ver/')
