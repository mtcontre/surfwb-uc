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
~~    * runs well with 1 processor~~
~~* save batifiles in output directory, and remove x,y,z from sol2d~~
~~* fix ver.py to read this~~
~~* reduce cfl-dt~~
~~* assign CB=-1 in interior boundaries (decomp2d)~~
~~* communicate through exchng2d~~
    *it fails with P=16,17,18,19
    *fixed!!!! the trick was on moving all blocks in mpisurf_BCS.f90 like
    
  DO i=3,Nx+2
    !Metricas
    xit(:,i,Ny+4)=xi(:,i-2,Ny-1)	!xiNy+4=xiNy-1
    xit(:,i,Ny+3)=xi(:,i-2,Ny)	!xiNy+3=xiNy

    etat(:,i,Ny+4)=eta(:,i-2,Ny-1)	!etaNy+4=etaNy-1
    etat(:,i,Ny+3)=eta(:,i-2,Ny)	!etaNy+3=etaNy

    zt(i,Ny+4)=zt(i,Ny+1)
    zt(i,Ny+3)=zt(i,Ny+2)
  END DO  
  
    above the "call exchange_2d" line (l.113)
* ~~fix verbose (master only)~~
* ~~distribute friction matrix~~
* ~~re-write outputgauges.f90~~
*~~ save executions time and number of processes~~
* add tests: thacker & pascua  and measure speedup
* Thacker: cartesian+variablebathymetry+wet-dry
    *setrun.py ready
    *figures dont match after 0.2T, with different processors
    *maybe some issue with the minimum height??
    *create editsub.py and the pbs file
    *spam!!!
~~* pascua:~~

~~* it fails for no reason :S see test6. I run the thacker case in ~/Downloads/SurfWB-UC_PAR/ and it looks great~~
~~* i should try to clean the other version.~~

#New version
The old version (the one I sent to Carlos) was OK, except for the bug that Marite found.
So I modified its makefiles, fixed some .f90 files and copied everything in the bitbucket repo.

#New TODO
on the Thacker case
* fix mpisurf_input_control.f90 so it matches the same input.dat as in the sequential version
* fixed timestep feature
* outputgauges.f90 in parallel
* measure new speed-up
* cant think of anything else -_-

