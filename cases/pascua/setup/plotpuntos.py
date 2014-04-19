from pylab import *
s=loadtxt('../data/bathy.dat')
x=reshape(s[:,0],(130,30))
y=reshape(s[:,1],(130,30))
z=reshape(s[:,2],(130,30))
b=loadtxt('../data/gauges.dat',skiprows=1)

f=figure(figsize=(16.,5))
pcolor(x,y,z,edgecolor='k')
scatter(b[:,1],b[:,2],alpha=0.5,linewidth=0.5,color='b')
axis('tight')
f.tight_layout()
savefig('puntos.png',dpi=600)