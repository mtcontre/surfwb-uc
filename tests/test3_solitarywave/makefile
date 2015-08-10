#need $(SURF) to be declared on .bashrc file

SURF_IN_DIR=$(CURDIR)/data/
SURF_OUT_DIR=$(CURDIR)/results/

SURF_LIB=$(SURF)/source
SURF_MAKE=$(SURF_LIB)/makefile


SURF_SOURCES = $(SURF_LIB)/MODULES.f90 \
	$(SURF_LIB)/main.f90 \
	$(SURF_LIB)/init.f90 \
	$(SURF_LIB)/input_control.f90 \
	$(SURF_LIB)/input_bound.f90 \
	$(SURF_LIB)/readIO.f90 \
	./readGA.f90 \
	$(SURF_LIB)/input_geom.f90 \
	$(SURF_LIB)/input_ic.f90 \
	$(SURF_LIB)/init_TS.f90 \
	$(SURF_LIB)/metrics.f90 \
	$(SURF_LIB)/Angulos.f90 \
	$(SURF_LIB)/massbalance.f90 \
	$(SURF_LIB)/tstep.f90 \
	$(SURF_LIB)/solver.f90 \
	$(SURF_LIB)/solver2.f90 \
	$(SURF_LIB)/solverf2.f90 \
	$(SURF_LIB)/solverf4.f90 \
	$(SURF_LIB)/friccion.f90 \
	$(SURF_LIB)/BCS.f90 \
	$(SURF_LIB)/fluxes.f90 \
	$(SURF_LIB)/vfroe1.f90 \
	$(SURF_LIB)/cfl_celerities.f90 \
	$(SURF_LIB)/outputmat.f90 \
	$(SURF_LIB)/outputgauges.f90 \
	$(SURF_LIB)/genabs0xi.f90 \
	$(SURF_LIB)/genabs0xi_2.f90 \
	$(SURF_LIB)/genabsNxi.f90 \
	$(SURF_LIB)/genabsNxi_2.f90 \
	$(SURF_LIB)/Inflow_xi.f90 \
	$(SURF_LIB)/Inflow_eta.f90 \
	$(SURF_LIB)/Outflow_xi.f90 \
 	$(SURF_LIB)/interp1.f90  \
	$(SURF_LIB)/tau.f90 \
 	$(SURF_LIB)/fzero_B.f90 \
 	$(SURF_LIB)/fzero_B2.f90  \
 	$(SURF_LIB)/fzeroEP_1.f90 \
 	$(SURF_LIB)/q0.f90 
	
include $(SURF_MAKE)

runsurf: xsurf
	./xsurf
	
cleanresults:
	rm -f $(SURF_OUT_DIR)*.gz
	rm -f $(SURF_OUT_DIR)*.dat