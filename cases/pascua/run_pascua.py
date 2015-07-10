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
    
def makegauges(indir):
  b=[[0, 1.00000, 1.00000],\
    [24, 2.17500, 0.72000 ],\
  [23, 2.51700, 0.79000 ],\
  [22, 3.00300, 0.81000 ],\
  [21, 3.61900, 0.78000 ],\
  [20, 4.22400, 0.83000 ],\
  [19, 4.83400, 0.95000 ],\
  [18, 5.66500, 1.04000 ],\
  [17, 6.27300, 1.03000 ],\
  [16, 6.88400, 1.10000 ],\
  [15, 7.49000, 1.50000 ],\
  [14, 8.10100, 1.92000 ],\
  [13, 8.71500, 2.04000 ],\
  [12, 9.33100, 2.05000 ],\
  [11, 9.93500, 2.21000 ],\
  [10, 10.57200, 2.42000 ],\
  [9, 11.17800, 2.68000 ],\
  [8, 11.79100, 2.66000 ],\
  [7, 12.40100, 2.44000 ],\
  [6, 13.01100, 2.38000 ]]
  f=open('%s/gauges.dat'%indir,'w')
  f.write('%i\n'%len(b))
  for i in range(len(b)):
    f.write('%i %.18e %.18e\n'%(int(b[i][0]),b[i][1],b[i][2]))
  f.close()

def testgauges(indir):
  g = np.loadtxt('%s/gauges.dat'%indir,skiprows=1)
  x = np.loadtxt('%s/gridx.dat'%indir)
  y = np.loadtxt('%s/gridy.dat'%indir)
  z = np.loadtxt('%s/gridz.dat'%indir)
  
  plt.figure(figsize=(9.,4.))
  plt.pcolormesh(x,y,z,cmap=plt.cm.gray,vmin=0.3,vmax=1.2)
  plt.colorbar()
  plt.axis('equal')
  plt.scatter(g[:,1],g[:,2],s=10.,facecolor='w',edgecolor='none')
  plt.xlim(x.min(),x.max())
  for i in range(g.shape[0]):
    plt.text(g[i,1],g[i,2],'%i'%g[i,0],fontsize=8,color='w')  
  plt.savefig('%s/gauges.png'%indir,bbox_inches='tight')  
  plt.close()
  
#--------parameters-------
def setrun(nxi,neta,nn,npr,indir):
  batifiles = ['%s/gridx.dat'%indir,'%s/gridy.dat'%indir,'%s/gridz.dat'%indir]
  initqfiles=['%s/inith.dat'%indir, '%s/initu.dat'%indir, '%s/initv.dat'%indir]
  caso=999
  tinit=0.
  tfinal=10.
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
  makegauges('data')
  testgauges('data')
  

