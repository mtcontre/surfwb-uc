#need $(SURF) to be declared on .bashrc file

F90 = mpif90
SURF_LIB=$(SURF)/source
SURF_MAKE=$(SURF_LIB)/makefile

SURF_COMMON = $(SURF_LIB)/input_bound.f90 \
	$(SURF_LIB)/metrics.f90 \
	$(SURF_LIB)/solver.f90 \
	$(SURF_LIB)/solver2.f90 \
	$(SURF_LIB)/solverf2.f90 \
	$(SURF_LIB)/solverf4.f90 \
	$(SURF_LIB)/friccion.f90 \
	$(SURF_LIB)/fluxes.f90 \
	$(SURF_LIB)/vfroe1.f90 \
	./readGA.f90 \
	$(SURF_LIB)/readIO.f90 \
	$(SURF_LIB)/genabs0xi_2.f90 \
	$(SURF_LIB)/genabs0xi.f90 \
	$(SURF_LIB)/genabsNxi_2.f90 \
	$(SURF_LIB)/genabsNxi.f90 \
	$(SURF_LIB)/interp1.f90 \
	$(SURF_LIB)/fzero_B.f90 \
	$(SURF_LIB)/fzero_B2.f90 \
	$(SURF_LIB)/fzeroEP_1.f90 \
	$(SURF_LIB)/tau.f90 \
	$(SURF_LIB)/Angulos.f90	

include $(SURF_MAKE)