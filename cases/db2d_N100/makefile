#compile Surf-UC
F90=mpif90
SURF_IN_DIR=$(CURDIR)/data/
SURF_OUT_DIR=$(CURDIR)/results/
SURF_PAR = /home/jgm/Downloads/SurfWB-UC_PAR/
SURF_LIB=$(SURF_PAR)/source/
SURF_MAKE=$(SURF_LIB)/makefile

SURF_SOURCES ?= $(SURF_LIB)MODULES.f90\
	$(SURF_LIB)main.f90 \
	$(SURF_LIB)init.f90 \
	$(SURF_LIB)input_control.f90 \
	$(SURF_LIB)input_bound.f90 \
	$(SURF_LIB)input_friction.f90 \
	$(SURF_LIB)readGA.f90 \
	$(SURF_LIB)readIO.f90 \
	$(SURF_LIB)input_geom.f90 \
	$(SURF_LIB)input_ic.f90 \
	$(SURF_LIB)init_TS.f90 \
	$(SURF_LIB)decomp_2d.f90 \
	$(SURF_LIB)bcast_bcs.f90 \
	$(SURF_LIB)metrics.f90 \
	$(SURF_LIB)Angulos.f90 \
	$(SURF_LIB)cfl_celerities.f90 \
	$(SURF_LIB)outputmat_par.f90 \
	$(SURF_LIB)verbs.f90 \
	$(SURF_LIB)massbalance.f90 \
	$(SURF_LIB)setdt.f90 \
	$(SURF_LIB)solver.f90 \
	$(SURF_LIB)solver2.f90 \
	$(SURF_LIB)solverf2.f90 \
	$(SURF_LIB)solverf4.f90 \
	$(SURF_LIB)friccion.f90 \
	$(SURF_LIB)BCS.f90 \
	$(SURF_LIB)exchange_2d.f90 \
	$(SURF_LIB)fluxes.f90 \
	$(SURF_LIB)genabs0xi_2.f90 \
	$(SURF_LIB)genabs0xi.f90 \
	$(SURF_LIB)genabsNxi_2.f90 \
	$(SURF_LIB)genabsNxi.f90 \
	$(SURF_LIB)Inflow_eta.f90 \
	$(SURF_LIB)Inflow_xi.f90 \
	$(SURF_LIB)Outflow_xi.f90 \
	$(SURF_LIB)interp1.f90\
	$(SURF_LIB)vfroe1.f90 \
	$(SURF_LIB)fzero_B.f90 \
	$(SURF_LIB)fzero_B2.f90 \
	$(SURF_LIB)fzeroEP_1.f90 \
	$(SURF_LIB)q0.f90 \
	$(SURF_LIB)tau.f90 

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
