#need $(SURF) to be declared on .bashrc file

F90 = mpif90
INDIR=$(CURDIR)/data
SURF_LIB=$(SURF)/source
SURF_MAKE=$(SURF_LIB)/makefile


include $(SURF_MAKE)


runsurf: xsurf
	./xsurf
	
cleanresults:
	rm -f $(SURF_OUT_DIR)*.gz
	rm -f $(SURF_OUT_DIR)*.dat
	
	