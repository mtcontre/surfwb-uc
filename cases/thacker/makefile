#compile Surf-UC
F90=mpif90
SURF_IN_DIR=$(CURDIR)/data/
SURF_OUT_DIR=$(CURDIR)/results/
SURF_PAR = /home/jgm/Downloads/SurfWB-UC_PAR/
SURF_LIB=$(SURF_PAR)/source/
SURF_MAKE=$(SURF_LIB)/makefile


include $(SURF_MAKE)

runsurf: xsurf
	./xsurf
	
cleanr:
	rm -f $(SURF_OUT_DIR)*.gz
	rm -f $(SURF_OUT_DIR)*.dat
run:	xsurf
	mpiexec.openmpi -np 4 xsurf
plots:
	rm -f visualizar/*.png
	cd visualizar;mpiexec.openmpi -np 3 python ver_par.py
	eog visualizar/frame*0.png
cleanplots:
	rm -f visualizar/*.png
cleanmazinger:
	rm -f *.o*
	rm -f *.e*
	rm -f 
	rm -f proj_nodes*
	rm -f -r N*
#
