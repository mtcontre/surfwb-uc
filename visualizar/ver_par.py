import pypar
from pylab import loadtxt,reshape,\
  pcolormesh,plot,show,colorbar,title,arange,savefig,close
from string import split
l=pypar.size()
r=pypar.rank()
print 'this is rank %i of %i '%(r,l)
# Read general parameters
f=open('../results/grids/gridproperties.dat')

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
nxi=line[1:]
nxi=[int(nxi[i]) for i in range(ngrids)]

line=f.readline();line=split(line)
neta=line[1:]
neta=[int(neta[i]) for i in range(ngrids)]

#set the range of frames to plot
t=loadtxt('../results/time.dat')
frames=range(0,len(t)*dit,dit)
rango=arange(0,len(frames),l)+r

batinames=[]
for lgrid in range(ngrids):
  line=f.readline();line=split(line)
  batinames.extend([line])
  batinames[lgrid]=[split(batinames[lgrid][i],sep='/')[1] for i in range(3)]
f.close()

for i in rango:
  if i<=len(frames)-1:
    it=frames[i]
  else:
    break
  print it
  
  #get bathymetry filenames  for each leve  
  for l in range(ngrids):
    x=loadtxt('../results/'+batinames[l][0])
    y=loadtxt('../results/'+batinames[l][1])
    z=loadtxt('../results/'+batinames[l][2])
    for i in range(dims[0]):
      for j in range(dims[1]):
	#file with properties for this grid
	fname='../results/grids/grid%03d_%03d.dat'%(i,j)
	f=open(fname)
	line=f.readline()
	nbx=int(split(line)[0]);line=f.readline()
	nby=int(split(line)[0]);line=f.readline()
	si=int(split(line)[0])-1;line=f.readline()#minus one for python
	ei=int(split(line)[0])+1;line=f.readline()
	sj=int(split(line)[0])-1;line=f.readline()
	ej=int(split(line)[0])+1;line=f.readline()
	#coord1=int(split(line)[0]);line=f.readline()
	#coord2=int(split(line)[0]);line=f.readline()
	f.close()
	fname='../results/frame%08d.%03d_%03d.dat'%(it,i,j)
	s=loadtxt(fname)
	h=reshape(s[:,0],(nbx,nby))
	pcolormesh(x[si:ei,sj:ej],y[si:ei,sj:ej],h,vmin=0.,vmax=10.)#x[si:ei,sj:ej],y[si:ei,sj:ej],
	ei=ei-2
	ej=ej-2
	#colorbar()
	plot([x[si,sj],x[si,ej],x[ei,ej],x[ei,sj],x[si,sj]],\
	  [y[si,sj],y[si,ej],y[ei,ej],y[ei,sj],y[si,sj]],color='k')

	    #linewidth=3,color='k')
	title('t=%.3fs'%(t[it/dit]))      
  savefig('frame%08d.png'%it)
  close()