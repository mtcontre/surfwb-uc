import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
T = 2.24
a = 1.0
r0 = 0.8
h0 = 0.1
g = 9.81
omega = np.sqrt(8*g*h0/a**2)
#A = (a**2-h0**2)/(a**2+h0**2)
A = (a**4-r0**4)/(a**4+r0**4)

def bati(x,y):
  r = np.sqrt(x**2+y**2)
  return -h0*(1-r**2/a**2)

def analytical(t,x,y):
  r = np.sqrt(x**2+y**2)
  z = -h0*(1-r**2/a**2)
  eta = h0*(np.sqrt(1-A**2)/(1-A*np.cos(omega*t))-1. \
	    -r**2/a**2*((1.-A**2)/(1.-A*np.cos(omega*t))**2-1.))
  u  = 1./(1.-A*np.cos(omega*t))*0.5*omega*x*A*np.sin(omega*t)
  v  = 1./(1.-A*np.cos(omega*t))*0.5*omega*y*A*np.sin(omega*t)
  
  #eta = np.where(eta<z, z, eta)
  #u = np.where(eta<z, 0., u)
  ##v = np.where(eta<z, 0., v)
  return eta,u,v
def eta0(x,y):
  eta,u,v = analytical(0.,x,y)
  return eta,u,v

def setrun(nxi,neta,nn,npr):
  #=====================================
  #input.dat parameters
  #=====================================
  import os  
  indir = 'data_%(nm)i_%(nn)i_%(np)i/'%({'nn':nn,'np':npr,'nm':nxi})
  outdir='results_%(nm)i_%(nn)i_%(np)i/'%({'nn':nn,'np':npr,'nm':nxi})
  os.system('mkdir %s'%indir)
  os.system('mkdir %s'%outdir)
  os.system('mkdir %s/timeseries'%outdir)  
  os.system('mkdir %s/grids/'%outdir)
  batifiles = ['%s/gridx.dat'%indir,'%s/gridy.dat'%indir,'%s/gridz.dat'%indir]
  initqfiles=['%s/inith.dat'%indir, '%s/initu.dat'%indir, '%s/initv.dat'%indir]
  
  caso = 999
  tinit = 0.0
  tfinal = T
  cfl = 0.7
  #nxi=x.shape[0]
  #neta=x.shape[1]
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
    dtout = T/5.
  kappa=1e-5
  rktype=1
  limtype=1
  fricopt=0
  outopt=1
  

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
    f.write('\n%.3f'%dtout)
  f.write('\n%3.20e'%kappa)
  f.write('\n%i'%rktype)
  f.write('\n%i'%limtype)
  f.write('\n%i'%fricopt)
  f.write('\n%i'%outopt)
  f.write('\n%s\n'%outdir)
  f.close()

  import numpy as np
  import matplotlib.pyplot as plt


  


  #=====================================
  #computational domain
  #=====================================
  x=np.linspace(-1.3,1.3,nxi)
  y=np.linspace(-1.3,1.3,neta)
  x,y=np.meshgrid(x,y)
  z = bati(x,y)
  h,u,v = eta0(x,y)
  h = np.where(h-z>0,h-z,0.)


  plt.figure()
  plt.pcolormesh(x ,y ,z )
  plt.colorbar()
  plt.savefig('%s/bati.png'%indir)

  plt.figure()
  plt.pcolormesh(x ,y ,h )
  plt.colorbar()
  plt.savefig('%s/inith.png'%indir)

  plt.figure()
  plt.pcolormesh(x ,y ,h +z)
  plt.colorbar()
  plt.savefig('%s/initeta.png'%indir)

  plt.figure()
  plt.pcolormesh(x ,y ,u )
  plt.savefig('%s/initu.png'%indir)
  plt.colorbar()

  plt.figure()
  plt.pcolormesh(x ,y ,v )
  plt.savefig('%s/initv.png'%indir)
  plt.colorbar()


  np.savetxt('%s/gridx.dat'%indir,x )
  np.savetxt('%s/gridy.dat'%indir,y )
  np.savetxt('%s/gridz.dat'%indir,z )
  np.savetxt('%s/inith.dat'%indir,h )
  np.savetxt('%s/initu.dat'%indir,u )
  np.savetxt('%s/initv.dat'%indir,v )

  #=====================================
  #wave gauges
  #=====================================
  #b=[[1,-25., 0.],
    #[2,0., 0.],
    #[3,25., 0.]]; 
  b = []
  f=open('%s/gauges.dat'%indir,'w')
  f.write('%i\n'%len(b))
  #for i in range(len(b)):
    #f.write('%i %.18e %.18e\n'%(b[i][0],b[i][1],b[i][2]))
  f.close()
  return indir
if __name__=='__main__':
  import sys
  arg = sys.argv[1:]
  nxi = int(arg[0])
  p = int(arg[1])
  nn = 1

  indir = setrun(nxi, nxi, nn, p)