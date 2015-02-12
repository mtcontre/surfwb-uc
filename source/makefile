#SURFWB-UC default makefile

F90?= mpif90
FFLAGS=-O3 -g
LFLAGS = 
SURF_LIB ?= ./

SURF_SOURCES ?= $(SURF_LIB)/module_mpisurf.f90 \
	$(SURF_LIB)/module_globals.f90 \
	$(SURF_LIB)/module_geometries.f90 \
	$(SURF_LIB)/module_jacobianos.f90 \
	$(SURF_LIB)/main_surfpar.f90 \
	 $(SURF_LIB)/init_par.f90 \
	 $(SURF_LIB)/input_control.f90 \
	 $(SURF_LIB)/init_params.f90 \
	 $(SURF_LIB)/input_geom.f90 \
	 $(SURF_LIB)/input_ic.f90 \
	 $(SURF_LIB)/input_friction.f90 \
	 $(SURF_LIB)/set_carttopo.f90 \
	 $(SURF_LIB)/decomp_domain.f90 \
	 $(SURF_LIB)/decomp_bcs.f90 \
	 $(SURF_LIB)/saveout_parameters.f90 \
	 $(SURF_LIB)/init_seq.f90 \
	 $(SURF_LIB)/adimension.f90 \
	 $(SURF_LIB)/metrics.f90 \
	 $(SURF_LIB)/massbalance.f90 \
	 $(SURF_LIB)/outputmat_par.f90 \
	 $(SURF_LIB)/setdt.f90 \
	 $(SURF_LIB)/solver.f90 \
	 $(SURF_LIB)/BCS.f90 \
	 $(SURF_LIB)/fluxes.f90 \
	 $(SURF_LIB)/cfl_celerities.f90 \
	 $(SURF_LIB)/vfroe1.f90 \
	 $(SURF_LIB)/exchange_2d.f90
	 

SURF_OBJECTS = $(subst .f90,.o, $(SURF_SOURCES) )

.PHONY: clean test
xsurf: $(SURF_OBJECTS)
	$(F90) $(LFLAGS) $(SURF_OBJECTS) -o $(CURDIR)/$@

%.o %.mod: %.f90
	$(F90) $(FFLAGS) -c $< -o $@ -J $(SURF_LIB)/
clean:
	rm -f $(SURF_LIB)/$(SURF_LIB)*.o
	rm -f $(SURF_LIB)/$(SURF_LIB)*.o
	rm -f $(SURF_LIB)/$(SURF_LIB)*.mod
	rm -f $(SURF_LIB)/$(SURF_LIB)*.mod
	rm -f xsurf
test: 	
	@echo $(SURF_SOURCES)
	@echo $(SURF_MODULES)
	@echo $(SURF_OBJECTS)
	@echo $(SURF_MODS)
	