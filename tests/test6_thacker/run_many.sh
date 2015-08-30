for i in $(seq 1 1 25)
do
  echo running with $i procs
  python setrun.py 100 $i
  export INDIR='data_100_1_'$i
  echo input dir $INDIR
  mpirun -np $i xsurf  
  python ver.py 100 $i
done