#need $(SURF) to be declared on .bashrc file

SURF_IN_DIR=$(CURDIR)/data/
SURF_OUT_DIR=$(CURDIR)/results/

SURF_LIB=$(SURF)/source
SURF_MAKE=$(SURF_LIB)/makefile


F90=mpif90

include $(SURF_MAKE)


runsurf: xsurf
	./xsurf
	
cleanresults:
	rm -f $(SURF_OUT_DIR)*.gz
	rm -f $(SURF_OUT_DIR)*.dat