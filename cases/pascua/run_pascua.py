import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt


def maketopo():
  s = np.loadtxt('SOL2D.0.dat')
  x = np.reshape(s[:,0],(130,30),order='F')
  y = np.reshape(s[:,1],(130,30),order='F')
  z = np.reshape(s[:,2],(130,30),order='F')
  h = np.reshape(s[:,3],(130,30),order='F')
  u = np.reshape(s[:,4],(130,30),order='F')
  v = np.reshape(s[:,5],(130,30),order='F')
  np.savetxt('data/gridx.dat',x)
  np.savetxt('data/gridy.dat',y)
  np.savetxt('data/gridz.dat',z)
  np.savetxt('data/inith.dat',h)
  np.savetxt('data/initu.dat',u)
  np.savetxt('data/initv.dat',v)
  

def test_topo():
    x = np.loadtxt('data/gridx.dat')
    y = np.loadtxt('data/gridy.dat')
    z = np.loadtxt('data/gridz.dat')
    h = np.loadtxt('data/inith.dat')
    u = np.loadtxt('data/initu.dat')
    v = np.loadtxt('data/initv.dat')
    
    plt.figure()
    plt.pcolormesh(x,y,z)
    plt.colorbar()
    plt.savefig('data/topo.png')
    plt.close()
    
    plt.figure()
    plt.pcolormesh(x,y,h)
    plt.colorbar()
    plt.savefig('data/inith.png')
    plt.close()
    
    plt.figure()
    plt.pcolormesh(x,y,u)
    plt.colorbar()
    plt.savefig('data/initu.png')
    plt.close()
    
    plt.figure()
    plt.pcolormesh(x,y,v)
    plt.colorbar()
    plt.savefig('data/initv.png')
    plt.close()
    

#--------parameters-------
def setrun(nxi,neta,nn,npr,indir):
  batifiles = ['%s/gridx.dat'%indir,'%s/gridy.dat'%indir,'%s/gridx.dat'%indir]
  initqfiles=['%s/inith.dat'%indir, '%s/initu.dat'%indir, '%s/initv.dat'%indir]
  caso=999
  tinit=0.
  tfinal=1.
  cfl=0.95
  #nxi=300
  #neta=300
  batiopt=1
  initqopt=1
  dxi=1.
  deta=1.
  L=1.
  H=1.
  U=1.
  bcxi0=1
  bcxiN=1
  bceta0=1
  bcetaN=1
  dit=-1
  if dit==-1:
    dtout=1.
  kappa=1e-5
  rktype=1
  limtype=1
  fricopt=0
  outopt=1
  outdir='results'#_%(nn)i_%(np)i_%(nm)i/'%({'nn':nn+1,'np':npr,'nm':nxi})

  import os
  os.system('mkdir %s'%outdir)
  os.system('mkdir %s/timeseries'%outdir)
  os.system('mkdir %s'%indir)
  os.system('export INDIR=%s'%indir)
  os.system('export OUTDIR=%s'%outdir)
  #---------write to file
  f=open('%s/input.dat'%indir,'w')
  f.write('%i'%caso)
  f.write('\n%.8f'%tinit)
  f.write('\n%.8f'%tfinal)
  f.write('\n%.8f'%cfl)
  f.write('\n%i'%nxi)
  f.write('\n%i'%neta)
  f.write('\n%i'%batiopt)
  for i in range(len(batifiles)):
    f.write('\n%s'%batifiles[i])
  f.write('\n%i'%initqopt)
  for i in range(len(initqfiles)):
    f.write('\n%s'%initqfiles[i])
  f.write('\n%.8f'%dxi)
  f.write('\n%.8f'%deta)
  f.write('\n%.8f'%L)
  f.write('\n%.8f'%H)
  f.write('\n%.8f'%U)
  f.write('\n%i'%bcxi0)
  if bcxi0==4:
    f.write('\n%i'%GA1)
    if GA1==9 or GA1==1:
      f.write('\n%i'%Nsenal1)
  f.write('\n%i'%bcxiN)
  f.write('\n%i'%bceta0)
  f.write('\n%i'%bcetaN)
  f.write('\n%i'%dit)
  if dit==-1:
    f.write('\n%.18e'%dtout)
  f.write('\n%.18e'%kappa)
  f.write('\n%i'%rktype)
  f.write('\n%i'%limtype)
  f.write('\n%i'%fricopt)
  f.write('\n%i'%outopt)
  f.write('\n%s\n'%outdir)
  f.close()
  print '---%s/input.dat generated---'%(indir)


if __name__=='__main__':
  setrun(130,30,1,1,'data')
  maketopo()
  test_topo()
  

