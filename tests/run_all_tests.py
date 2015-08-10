import os
os.chdir('test1_db1d')
os.system('rm -r data results vis')
os.system('rm xsurf')
os.system('make')
os.system('python setrun.py')
os.environ['INDIR'] = 'data/'
os.system('echo $INDIR')
os.system('./xsurf')