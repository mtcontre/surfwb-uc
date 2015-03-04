def writesub(nn,np,nm, memory, indir):
  """
    nn = number of nodes
    np = number of processors per node
    nm = number of mesh divisions in the anuga.create_rectangular_cross_domain
  """
  nn=int(nn)
  np=int(np)
  s="""#!/bin/bash
#PBS -N db2d_%(nn)i_%(np)i_%(nm)i
#PBS -l nodes=%(nn)i:ppn=%(np)i
#PBS -l walltime=01:59:59
#PBS -l pmem=%(mb)iMB
#PBS -M jdgalaz@uc.cl
#PBS -V
#PBS -q work
#PBS -m ea
cd $PBS_O_WORKDIR
date
ulimit -s unlimited
module load openmpi/1.6.5/gcc
module load anaconda/1.8.0
export OMP_STACKSIZE=16M
export OMP_NUM_THREADS=%(np)i
make
export INDIR=%(indir)s
./xsurf
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
  from run_db2d import setrun

  ies =[[1,4,8,16,32]]
  nmesh = [100,400,600,1000]
  memory = [1000,1000,1000,2000.]
  nnodes = 1

  #dividir el loop: n para nodos, i= proc por nodo
  for nm,mem in zip(nmesh,memory):
   for n in range(nnodes):
    for i in ies[n]:     
      #create input.dat and bathy-initq
      indir='data_%(nn)i_%(np)i_%(nm)i/'%({'nn':n+1,'np':i,'nm':nm})
      setrun(nm, nm, n, i, indir)
      
      #write pbs file and send to queue
      mb = min([mem, int(32000./i)])
      writesub(n+1,i,nm, mb, indir)      
      #msg=os.system('qsub nsub.sub')
      #time.sleep(2)

