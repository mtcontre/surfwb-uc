for i in $(seq 1 1 10)
do
  mpirun -np $i xsurf
  python ver.py
done