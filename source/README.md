Structure of source directory:

* oldsources contains old .f90 files that must be checked for deletion
* seqsource contains .f90 sources that run exclusively by the sequential form of the code (main_seq.f90 genabs0xi.f90 for instance)
* parsource: the same but mpiparallel (exchange2d.f90 for example)
* $SURF/source contains common routines to both implementations