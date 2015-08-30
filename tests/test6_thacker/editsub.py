def writesub(nn,np,nm, memory, indir):
  """
    nn = number of nodes
    np = number of processors per node
    nm = number of mesh divisions in the anuga.create_rectangular_cross_domain
  """
  nn=int(nn)
  np=int(np)
  s="""#!/bin/bash
#PBS -N thk_%(nm)i_%(nn)i_%(np)i
#PBS -l nodes=%(nn)i:ppn=%(np)i
#PBS -l walltime=01:59:59
#PBS -l pmem=%(mb)iMB
#PBS -M jdgalaz@uc.cl
#PBS -V
#PBS -q test
#PBS -m ea
cd $PBS_O_WORKDIR
date
ulimit -s unlimited
module load openmpi/1.6.5/gcc
module load anaconda/1.8.0
make
python setrun.py %(nm)i %(np)i
export INDIR=%(indir)s
mpirun xsurf
python ver.py
date"""%({'nn': nn, 'np': np, 'nm': nm, 'mb': memory, 'indir':indir})
  f=open('nsub.sub','w')
  f.write(s)
  f.close()
def qdel(n0,nf):
  import os
  for n in range(n0,nf+1):
    os.system('qdel %i'%n)

if __name__ == '__main__':
  import os
  import time
  from setrun import setrun

  pes =[[1,2,4,8,16,32]]
  nmesh = [100]
  memory = [1000]
  nnodes = 1
  #dividir el loop: n para nodos, i= proc por nodo
  for nm,mem in zip(nmesh,memory):
   for n in range(nnodes):
    for p in pes[n]:     
      #create input.dat and bathy-initq
      indir='data_%(nm)i_%(nn)i_%(np)i/'%({'nn':n+1,'np':p,'nm':nm})
      setrun(nm, nm, n+1, p, indir)
      
      #write pbs file
      mb = min([mem, int(32000./p)])
      writesub(n+1,p,nm, mb, indir)    
      
      #send to q
      os.system('qsub nsub.sub')
      
 
