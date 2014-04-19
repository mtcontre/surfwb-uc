#compile Surf-UC
F90=mpif90
SURF_IN_DIR=$(CURDIR)/data/
SURF_OUT_DIR=$(CURDIR)/results/

SURF_LIB=../../source/
SURF_MAKE=$(SURF_LIB)/makefile

include $(SURF_MAKE)

.PHONY: xsurf cleanresults cleanplot plot
run: xsurf
	./xsurf
run3: xsurf
	mpiexec.openmpi -np 3 xsurf
cleanresults:
	rm -f $(SURF_OUT_DIR)*.gz
	rm -f $(SURF_OUT_DIR)*.dat
	rm -f visualizar/*.png
cleanplot:
	rm visualizar/*.png
plot:
	cd visualizar
	mpiexec.openmpi -np 3 python ver_par.py
	