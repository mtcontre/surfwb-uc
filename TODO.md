#TODO list

##General fixes
* fix readGA.f90 to read h01 more smart
* check and remove sources/oldsources/*.f90

##Parallel dev
Steps:
Working on test2_db2d:

* ~~ separate makefile source lists in: parallel, sequential and common ~~
* ~~ add mpisurf_* on duplicate files ~~
* ~~ leave only periodic/open/closed boundary conditions in BCS.f90~~
* ~~ add init and finalize clauses to main.f90 ~~
* ~~ run code with 4 processes, each one saving files called SOL2D.P.it.dat.gz~~
* ~~Create cartesian topology~~
* ~~calculate si,ei,sj,ej, nbx,nby for each process~~
* ~~ determine which matrices should be extracted and red~~
    * ~~all arrays allocated in input_geom.f90, input_ic.f90~~
    * ~~no others array are declared before decomp2d in mpisurf_init.f90~~
  
~~* modify outputmat and outputgauges to generate parallel output properly~~
~~    *mpisurf_outputgauges modification will be after master merging,now commented from main!!~~
~~    *mpisurf_outputmat was ready since before! :)~~
~~* save topology information to a file for later check~~
~~    * save total nxi and neta (nbx,nby before decomp2d) in variable old_nbx,old_nby~~
~~* run test2 with changes up to this point~~
    * runs well with 1 processor
* save batifiles in output directory, and remove x,y,z from sol2d
* fix ver.py to read this
* reduce dt
* communicate