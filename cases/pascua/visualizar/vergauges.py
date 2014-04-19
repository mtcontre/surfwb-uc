def gaugessurfuc2(fname):
#lee el 'fname', que es un fort.gauge salido del geoclaw
#retorna los id de los gauges =gaugenos
#retorna el vector tiempo para los gauges = time_clawcart
#retorna la matriz de los eta(g,t) =eta_clawcart
  from matplotlib.mlab import find
  def reordenar(eta,orden):
    eta1=1*eta
    for i in range(len(orden)):
      eta1[i]=eta[orden[i]]
    return eta1

  def stars2num(s):
      if s[0]=='*':
	gaugeno=99999
      else:
	gaugeno=int(s)
      return gaugeno  
  gdata=np.loadtxt(fname,converters={0:stars2num})  
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
  orden=[0,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
  print len(gaugenos)
  print orden
  gaugenos=reordenar(list(gaugenos),orden)
  eta_clawcart=reordenar(eta_clawcart,orden)
  return gaugenos,eta_clawcart,time_clawcart