
include ../site_specific.mk_${HOSTTYPE}

default: without_mpi

modules:
	cd etc; make timing_module.o
	cd density; make density.o
	cd main; make main_module.o
	cd idemix; make idemix_module.o
	cd diagnostics; make diagnostics_module.o

all:    modules
	cd main; make 
	cd idemix; make
	cd density; make
	cd diagnostics; make 

config.o: all config.f90 
	$(F90) $(F90FLAGS) $(CDFFLAGS) -Imain  -Iidemix -Idiagnostics -Idensity -c config.f90

without_mpi: all config.o
	cd parallel; make parallel_none.o
	$(F90) main/*.o idemix/*.o  density/*.o diagnostics/*.o \
	  config.o parallel/parallel_none.o etc/timing_module.o\
	 $(F90FLAGS) $(CDFFLAGS)-o ../bin/model.x

dirs = density diagnostics etc idemix  main parallel 
clean: 
	for d in $(dirs); do cd $$d && make clean && cd ..; done
	rm -f *.o *.mod 


