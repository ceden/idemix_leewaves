
include ../../site_specific.mk_${HOSTTYPE}

all: main_module.o numerics.o main.o \

clean :
	 rm -f *.o *.mod 

main_module.o : main_module.f90 
	$(F90) $(F90FLAGS) -c main_module.f90
numerics.o: numerics.f90 main_module.o 
	$(F90) $(F90FLAGS) -c numerics.f90
main.o: main.f90 main_module.o ../idemix/idemix_module.o ../etc/timing_module.o
	$(F90) $(F90FLAGS) -I../idemix -I../diagnostics -I../etc -c main.f90


