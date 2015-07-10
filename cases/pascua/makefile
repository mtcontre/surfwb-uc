#compile Surf-UC
F90=gfortran
SURF_IN_DIR=$(CURDIR)/data/
SURF_OUT_DIR=$(CURDIR)/results/

SURF_LIB=../../source/
SURF_MAKE=$(SURF_LIB)/makefile
	
include $(SURF_MAKE)

.PHONY: xsurf cleanresults cleanplot plot
run: xsurf
	./xsurf
cleanloc:
	rm -f $(INDIR)
	rm -f results
	rm -f vis	