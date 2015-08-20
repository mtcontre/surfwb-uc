def writesub(nn,np):
  nn=int(nn)
  np=int(np)
  s=''
  s+='#!/bin/bash'
  s+='\n'+'#PBS -N db2d'
  s+='\n'+'#PBS -l nodes=%i:ppn=%i'%(nn,np)
  s+='\n'+'#PBS -l walltime=00:59:59'
  s+='\n'+'#PBS -l pmem=200MB'
  s+='\n'+'#PBS -M jdgalaz@uc.cl'
  s+='\n'+'#PBS -V'
  s+='\n'+'#PBS -q work'
  s+='\n'+'cd $PBS_O_WORKDIR'
  s+='\n'+'date'
  s+='\n'+'# For future check to record the names of nodes used for your job'
  s+='\n'+'sort -u -t- -n --key=2 --key=3 -u $PBS_NODEFILE > proj_nodes.$PBS_JOBID'
  s+='\n'+'ulimit -s unlimited'
  s+='\n'+'mpiexec xsurf'
  s+='\n'+'date'

  f=open('nsub.sub','w')
  f.write(s)
  f.close()
  
if __name__ == '__main__':
  import sys
  import os
  
  ies=[[1,8,16,24,32],[20,24,28,32],[24,27,30,32]]
  for n in range(3):
   for i in ies[n]:
    writesub(n+1,i)
    msg=os.system('qsub nsub.sub')
    print msg, (n+1)*i,n,i
    #f=open('log','a')
    #f.write('%s,%s,%s\n'%(msg,sys.argv[1],sys.argv[2]))
   #print sys.argv
