import numpy as np
import matplotlib.pyplot as plt
t_par1 = np.loadtxt('times.txt',delimiter=',')
t_par2 = np.loadtxt('times_omp2.txt',delimiter=',')

grids = [100., 400., 600., 1000.]

f1 = plt.figure(figsize=(16.,6.))

ax1 = f1.add_subplot(121)
ax2 = f1.add_subplot(122)

ccycle = plt.rcParams['axes.color_cycle']
for g,ig in zip(grids,range(len(grids))):
  i1 = np.where(t_par1[:,0]==g)[0]
  i2 = np.where(t_par2[:,0]==g)[0]
  
  p1 = t_par1[i1,1]
  p2 = t_par2[i2,1]
  
  #calc times
  t1 = t_par1[i1,2]
  t2 = t_par2[i2,2]
  
  ax1.plot(p1,t1,'o',label='%i'%g, color=ccycle[ig])
  ax1.plot(p2,t2,'s',label='%i'%g, color=ccycle[ig])
  
  ax1.set_yscale('log')   
  ax1.set_title('T')
  ax1.legend(loc=0,ncol=2)
  ax1.set_ylim(0.5e-1,1e5)
  #ax2.legend()
  
  #speedups
  #index of T1 for each case
  it1one = np.where(p1==1.)[0]
  it2one = np.where(p2==1.)[0]
  if len(it1one)>0:
    tone1 = t1[it1one]
    sp1 = tone1/t1
    ax2.plot(p1,sp1,'o',label='%i'%g, color=ccycle[ig])
  else:
    tone1 = None

  if len(it2one)>0:
    tone2 = t2[it2one]
    sp2 = tone2/t2
    print g,p2,sp2
    ax2.plot(p2,sp2,'s',label='%i'%g, color=ccycle[ig])
  else:
    tone2 = None
  ax2.set_title('Sp')
  
  
  
plt.show()