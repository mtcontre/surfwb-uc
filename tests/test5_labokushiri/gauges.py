import numpy as np
import matplotlib.pyplot as plt

def gaugessurfuc2(fname):
#lee el 'fname', que es un fort.gauge salido del geoclaw
#retorna los id de los gauges =gaugenos
#retorna el vector tiempo para los gauges = time_clawcart
#retorna la matriz de los eta(g,t) =eta_clawcart
  from matplotlib.mlab import find
  gdata=np.loadtxt(fname)  
  gaugeno=np.array(gdata[:,0],dtype=int)
  level=np.array(gdata[:,1],dtype=int)
  t=[]#np.zeros(gdata[:,2].shape)
  t=gdata[:,1]
  q=gdata[:,2:]
  gaugenos=set(gaugeno)
  ngauges=len(gaugenos)
  eta_clawcart=np.zeros((ngauges,len(t)/ngauges))
  time_clawcart=np.zeros(len(t)/ngauges)
  for n in range(len(gaugenos)):
    #n=int(n)
    nn = find(gaugeno==list(gaugenos)[n]) # vector con los indices correspondientes
    #gauges[n]=GaugeSolution()
    eta_clawcart[n]=q[nn,3]
    #print q[nn,:]
    if n==0:
      time_clawcart=t[nn]
  #orden=[0,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
  #print len(gaugenos)
  #print orden
  #gaugenos=reordenar(list(gaugenos),orden)
  #eta_clawcart=reordenar(eta_clawcart,orden)
  return gaugenos,eta_clawcart,time_clawcart

gauges,eta,time = gaugessurfuc2('results/timeseries/gauges.999.dat')
gauges=list(gauges)
f2=plt.figure()
#ax
for i in range(3):
  for j in range(3):
    k=3*(i)+j
    ax2=f2.add_subplot(3,3,k+1)
    ax2.plot(time,eta[k],color='r')
    ax2.set_title(gauges[k])
    ax2.set_xlim(0.,50.)