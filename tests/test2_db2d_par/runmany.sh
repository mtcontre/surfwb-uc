for i in $(seq 1 1 25)
do
  mpirun -np $i xsurf
  python ver.py
done