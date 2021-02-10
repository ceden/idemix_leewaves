
include ../../site_specific.mk_${HOSTTYPE}

all: diagnostics_module.o diag_main.o diag_snap_idemix.o 

clean :
	 rm -f *.o *.mod 

diagnostics_module.o: diagnostics_module.f90
	$(F90) $(F90FLAGS) -c diagnostics_module.f90

diag_main.o: diag_main.f90 ../main/main_module.o  ../idemix/idemix_module.o  diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../idemix  -c diag_main.f90


diag_snap_idemix.o: diag_snap_idemix.f90 ../main/main_module.o ../idemix/idemix_module.o diagnostics_module.o
	$(F90) $(F90FLAGS) $(CDFFLAGS) -I../main  -I../idemix  -c diag_snap_idemix.f90


