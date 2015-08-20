#need $(SURF) to be declared on .bashrc file
F90 = gfortran
INDIR=$(CURDIR)/data
SURF_LIB=$(SURF)/source
SURF_MAKE=$(SURF_LIB)/makefile

SURF_SEQ ?= $(SURF_LIB)/main.f90 \
	$(SURF_LIB)/init.f90 \
	$(SURF_LIB)/input_control.f90 \
	$(SURF_LIB)/input_geom.f90 \
	$(SURF_LIB)/input_ic.f90 \
	$(SURF_LIB)/input_bound.f90 \
	$(SURF_LIB)/BCS.f90 \
 	$(SURF_LIB)/outputmat.f90 \
	$(SURF_LIB)/outputgauges.f90 \
	$(SURF_LIB)/tstep.f90 \
	$(SURF_LIB)/Inflow_eta.f90 \
	$(SURF_LIB)/Inflow_xi.f90 \
	$(SURF_LIB)/Outflow_xi.f90 \
	$(SURF_LIB)/q0.f90 \
	$(SURF_LIB)/massbalance.f90 \
	$(SURF_LIB)/cfl_celerities.f90	

SURF_COMMON ?= $(SURF_LIB)/init_TS.f90 \
	$(SURF_LIB)/metrics.f90 \
	$(SURF_LIB)/solver.f90 \
	$(SURF_LIB)/solver2.f90 \
	$(SURF_LIB)/solverf2.f90 \
	$(SURF_LIB)/solverf4.f90 \
	$(SURF_LIB)/friccion.f90 \
	$(SURF_LIB)/fluxes.f90 \
	$(SURF_LIB)/vfroe1.f90 \
	$(SURF_LIB)/readGA.f90 \
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
	
SURF_PAR ?= $(SURF_LIB)/mpisurf_input_geom.f90 \
	$(SURF_LIB)/mpisurf_input_ic.f90 \
	$(SURF_LIB)/mpisurf_main.f90 \
	$(SURF_LIB)/mpisurf_init.f90 \
	$(SURF_LIB)/mpisurf_input_control.f90 \
	$(SURF_LIB)/mpisurf_input_friction.f90 \
	$(SURF_LIB)/mpisurf_decomp_2d.f90 \
	$(SURF_LIB)/mpisurf_bcast_bcs.f90 \
	$(SURF_LIB)/mpisurf_cfl_celerities.f90 \
	$(SURF_LIB)/mpisurf_outputmat.f90 \
	$(SURF_LIB)/mpisurf_verbs.f90 \
	$(SURF_LIB)/mpisurf_massbalance.f90 \
	$(SURF_LIB)/mpisurf_setdt.f90 \
	$(SURF_LIB)/mpisurf_BCS.f90 \
	$(SURF_LIB)/mpisurf_exchange_2d.f90
	
include $(SURF_MAKE)

run: xsurf
	./xsurf
	
cleanresults:
	rm -f $(SURF_OUT_DIR)*.gz
	rm -f $(SURF_OUT_DIR)*.dat