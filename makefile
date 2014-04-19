#Archivo para compilar y crear ejecutable de TMGP
# es comentario

F90	=mpif90 
#/opt/intel/Compiler/11.1/072/bin/intel64/ifort
SRC	=src
FFLAGS	=-O3 
LIBs	=/home/jose/tec360/lib/libtecio.a -lstdc++ #/opt/tecplot/lib/tecio64.a -lstdc++
# .PHONY = clean

.SUFFIXES: .f90 .o 
OBJECTS=  MODULES.o main.o init.o input_control.o input_bound.o input_friction.o \
	readGA.o readIO.o input_geom.o  input_ic.o init_TS.o decomp_2d.o bcast_bcs.o metrics.o Angulos.o \
	cfl_celerities.o outputmat_par.o verbs.o massbalance.o setdt.o solver.o \
	solver2.o solverf2.o solverf4.o friccion.o BCS.o exchange_2d.o fluxes.o genabs0xi_2.o genabs0xi.o \
	genabsNxi_2.o genabsNxi.o Inflow_eta.o Inflow_xi.o Outflow_xi.o interp1.o\
	vfroe1.o fzero_B.o fzero_B2.o fzeroEP_1.o q0.o tau.o 

#input_ic.o
xsurf: $(OBJECTS) makefile
#	$(F90) $(FFLAGS) -o $@ $(objects) $(LIBs)
	$(F90) $(FFLAGS) -o $@ $(OBJECTS)
	


.f90.o:
	$(F90)  $(FFLAGS) -c $*.f90

clean:
	rm -f $(OBJECTS) xsurf

	