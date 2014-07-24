#$ export ETS_TOOLKIT=wx
#$ ipython --gui=wx
from pylab import *
from matplotlib.colors import LinearSegmentedColormap, ColorConverter
black = [0.0,0.0,0.0]
white = [1.0,1.0,1.0]
red = [1.0,0.0,0.0]
green = [0.0,1.0,0.0];
dark_green = [0.1,0.4,0.0];
light_green = [0.8,1.0,0.5];
blue = [0.0,0.0,1.0];
dark_blue = [0.2,0.2,0.7];
light_blue = [0.5,0.5,1.0];
blue_green = [0.0,1.0,1.0];
tan = [0.9,0.8,0.2];
tan = [0.8,0.5,0.2];
brown = [0.9,0.8,0.2];
gray8 = [0.8,0.8,0.8];
purple = [0.8,0.3,0.8];
colors = {0.0:blue,
                                            0.5:blue_green,
                                            1.0:red}
def make_colormap(colors):                                         
  z = sort(colors.keys())                                            
  n = len(z)
  z1 = min(z)
  zn = max(z)
  x0 = (z - z1) / (zn - z1) 
  R = []
  G = []
  B = []
  for i in range(n):
      #i'th color at level z[i]:
      Ci = colors[z[i]]      
      if type(Ci) == str:
	  # a hex string of form '#ff0000' for example (for red)
	  RGB = CC.to_rgb(Ci)
      else:
	  # assume it's an RGB triple already:
	  RGB = Ci
      R.append(RGB[0])
      G.append(RGB[1])
      B.append(RGB[2])
      ind=int(x0[i])      
  cmap_dict = {}
  cmap_dict['red'] = [(x0[i],R[i],R[i]) for i in range(len(R))]
  cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
  cmap_dict['blue'] = [(x0[i],B[i],B[i]) for i in range(len(B))]
  mymap = LinearSegmentedColormap('mymap',cmap_dict)
  geo_water=mymap(arange(256))
  R=array([int(255*c) for c in geo_water.transpose()[0]])
  G=array([int(255*c) for c in geo_water.transpose()[1]])
  B=array([int(255*c) for c in geo_water.transpose()[2]])
  return R,G,B

geo_water=make_colormap(colors)
land1_colormap = make_colormap({0.0:dark_green,
                                          1000.0:green,
                                          2000.0:light_green,
                                          4000.0:tan})
                                         
land2_colormap = make_colormap({0:dark_green,
                                          50:green,
                                          100:light_green,
                                          200:tan})
                                          
water_land_colormap = make_colormap({-1000:dark_blue,
                                               -500:blue,
                                               0:light_blue,
                                               .1:tan,
                                               5:tan,
                                               6:dark_green,
                                               1000:green,
                                               2000:light_green,
                                               4000:tan})
                                               
bathy1_colormap = make_colormap({-1000:brown,
                                           0:tan,
                                           .1:dark_green,
                                           1000:green,
                                           2000:light_green})
       
bathy2_colormap = make_colormap({-1000:brown,
                                           -100:tan,
                                           0:dark_green,
                                           .1:dark_green,
                                           1000:green,
                                           2000:light_green})

bathy3_colormap = make_colormap({-1:[0.3,0.2,0.1],
                                           -0.01:[0.95,0.9,0.7],
                                           .01:[.5,.7,0],
                                           1:[.2,.5,.2]})

seafloor_colormap = make_colormap({-1:[0.3,0.2,0.1],
                                              0:[0.95,0.9,0.7]})

land_colormap = make_colormap({ 0:[0.95,0.9,0.7],
                                          1:[.2,.5,.2]})

land_colors = make_colormap({0:[.5,.7,0], 1:[.2,.5,.2]})
def test_colormap(cmap):
  from mayavi import mlab
  x,y=meshgrid(linspace(0,100,101),linspace(0,100,101))
  z=(x+y)*0.5
  if True:
      mesh_water=mlab.mesh(x.transpose(),y.transpose(),z.transpose(),colormap='gist_earth',vmin=0.,vmax=100.)
      #Retrieve the LUT of the surf object.
      lut = mesh_water.module_manager.scalar_lut_manager.lut.table.to_array()
      lut[:, 0] = cmap[0]
      lut[:, 1] = cmap[1]
      lut[:, 2] = cmap[2]
      mesh_water.module_manager.scalar_lut_manager.lut.table = lut
      mlab.draw()
      mlab.view(50.495536875966657,
 26.902697031665959,
 334.60652149512265,
 array([ 50.,  50.,  50.]))

      mlab.colorbar()
      
def test_all():
  from mayavi import mlab
  test_colormap(land_colormap)
  mlab.savefig('land_colormap.png')
  mlab.close()
  
  test_colormap(land_colors)
  mlab.savefig('land_colors.png')
  mlab.close()

  test_colormap(seafloor_colormap)
  mlab.savefig('seafloor_colormap.png')
  mlab.close()

  test_colormap(bathy3_colormap)
  mlab.savefig('bathy3_colormap.png')
  mlab.close()

  test_colormap(bathy2_colormap)
  mlab.savefig('bathy2_colormap.png')
  mlab.close()

  test_colormap(bathy1_colormap)
  mlab.savefig('bathy1_colormap.png')
  mlab.close()

  test_colormap(water_land_colormap)
  mlab.savefig('water_land_colormap.png')
  mlab.close()
  
  test_colormap(land2_colormap)
  mlab.savefig('land2_colormap.png')
  mlab.close()
  
  test_colormap(land1_colormap)
  mlab.savefig('land1_colormap.png')
  mlab.close()
  
  test_colormap(geo_water)
  mlab.savefig('geo_water.png')
  mlab.close()
  
