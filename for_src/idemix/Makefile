
include ../../site_specific.mk_${HOSTTYPE}

all: idemix_module.o idemix.o idemix_leewaves.o

clean :
	 rm -f *.o *.mod 

idemix_module.o  : idemix_module.f90 ../main/main_module.o 
	 	$(F90) $(F90FLAGS) -I../main -c idemix_module.f90
idemix.o : idemix.f90 ../main/main_module.o idemix_module.o
		$(F90) $(F90FLAGS) -I../main -c idemix.f90
idemix_leewaves.o : idemix_leewaves.f90 ../main/main_module.o idemix_module.o
		$(F90) $(F90FLAGS) -I../main -c idemix_leewaves.f90
