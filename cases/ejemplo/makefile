#compile Surf-UC
F90=mpif90
SURF_IN_DIR=$(CURDIR)/data/
SURF_OUT_DIR=$(CURDIR)/results/

SURF_LIB=../../source/
SURF_MAKE=$(SURF_LIB)/makefile

include $(SURF_MAKE)

runsurf: xsurf
	./xsurf
	
cleanresults:
	rm -f $(SURF_OUT_DIR)*.gz
	rm -f $(SURF_OUT_DIR)*.dat