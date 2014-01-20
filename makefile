#Archivo para compilar y crear ejecutable de TMGP
# es comentario

F90	=gfortran#/opt/intel/composer_xe_2013/bin/ifort 
#/opt/intel/Compiler/11.1/072/bin/intel64/ifort
SRC	=src
FFLAGS	=-O3 
LIBs	=/home/jose/tec360/lib/libtecio.a -lstdc++ #/opt/tecplot/lib/tecio64.a -lstdc++


.SUFFIXES: .f90 .o 
objects=  MODULES.o main.o init.o input_control.o input_bound.o readGA.o readIO.o \
	input_geom.o  input_ic.o init_TS.o metrics.o Angulos.o cfl_celerities.o \
	outputmat.o outputgauges.o massbalance.o tstep.o solver.o solver2.o \
	solverf2.o solverf4.o friccion.o BCS.o fluxes.o genabs0xi_2.o genabs0xi.o \
	genabsNxi_2.o genabsNxi.o Inflow_eta.o Inflow_xi.o Outflow_xi.o interp1.o\
	vfroe1.o fzero_B.o fzero_B2.o fzeroEP_1.o q0.o tau.o
# objects = MODULES.o main.o readGA.o readIO.o init.o input_control.o input_geom.o bati.o metrics.o \
# 	  BCS.o fluxes.o vfroe1.o solver.o outputmat.o outputgauges.o  massbalance.o \
# 	  solver2.o metrics8.o pascua.o solverf4.o friccion.o \
# 	  solverf2.o genabs0xi.o genabsNxi.o fzero.o interp1.o tau.o \
# 	  genabs0xi_2.o genabsNxi_2.o fzero_2.o PropSecNat.o \
# 	  Outflow_xi.o Inflow_xi.o Inflow_eta.o fzeroEP_1.o grado3.o Angulos.o CarSaliente.o \
# 	  fzero_B.o fzero_B2.o q0.o blend.o input_ic.o\
# 	  kdtree2.o LPT.o lpt_p.o lpt_c.o tracking.o locate.o init_TS.o input_bound.o cfl_celerities.o tstep.o
# 	  #genabs0eta.o genabsNeta.o 

xsurf: $(objects) makefile
#	$(F90) $(FFLAGS) -o $@ $(objects) $(LIBs)
	$(F90) $(FFLAGS) -o $@ $(objects)
	


.f90.o:
	$(F90) $(FFLAGS) -c $*.f90

cleanobj:
	rm *.o