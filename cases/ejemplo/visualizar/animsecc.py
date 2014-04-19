from pylab import *
from mayavi import mlab
#mlab.options.offscreen = True
#ion()
p=loadtxt('../results/param.dat')
case=int(p[0])
nx=int(p[1])
ny=int(p[2])
dit=int(p[-1])
t=loadtxt('../results/Time'+str(case)+'.dat')
ntoplot=len(t)-1#min(len(t)-1,10)
rango=linspace(0,ntoplot*dit,ntoplot+1)

x=array([])
y=array([])
#figure()
nfr=0
#pause
def analitic(x,t):
  #xvector
  #t scalar
  h=zeros(x.shape)
  u=zeros(x.shape)
  #print x.shape[1]
  for i in range(len(x)):
    if x[i]/t<-sqrt(9.81*10):
      h[i]=10
      u[i]=0
    elif x[i]/t <=2*sqrt(9.81*10):
      h[i]=1/(9*9.81)*(2*sqrt(9.81*10)-x[i]/t)**2
      u[i]=2/3*(x[i]/t+sqrt(9.81*10))
    else:
      h[i]=0
      u[i]=0
  return h,u
for i in rango:
  #i=100
  #print i
  nfr+=1
  i=int(i)
  s=loadtxt('../results/SOL2D.%i.dat' %i)
  print '%i'%i
  x=reshape(s[:,0],(ny,nx))
  y=reshape(s[:,1],(ny,nx))
  z=reshape(s[:,2],(ny,nx))
  h=reshape(s[:,3],(ny,nx))
  u=reshape(s[:,4],(ny,nx))
  v=reshape(s[:,5],(ny,nx))
  plot(x[5,:],h[5,:],label='calculado',marker='.')
  ha,ua=analitic(x[5,:],t[i/dit])
  plot(x[5,:],ha,label='sol.analitica')
  axis('tight')
  ylim(0,11) 
  legend(loc=0)
  title('T=%.4fs de %.4fs'%(t[i/dit],t[-1]))
  
  #show()
  savefig('sec'+str(i).rjust(6,'0')+'.png')
  close()
  
  #close()
  



# View it.

#s = 
