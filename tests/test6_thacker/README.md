I used this test mainly to measure parallel performance. 

To run a 100x100 grid with 4 processors use
    python setrun.py 100 4
    mpirun -np 4 xsurf

And to visualize results use
    python ver.py 100 4

this will create files in a directory called ver_100/