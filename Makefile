#  PGAS-FMO proxy app Makefile 

 FC=ftn         # theta
#FC=mpif90      # cooley
#FC=mpif90      # summit
 FLAGS=-O3 -w

all: pgas-fmo.x  

#pgas-fmo.f90: pgas-fmo.0.0.f90  
#	cp pgas-fmo.0.0.f90  pgas-fmo.f90 

pgas-fmo.x: src/pgas-fmo.f90  
	$(FC) $(FLAGS) -o bin/pgas-fmo.x src/pgas-fmo.f90 

clean:
	rm *.o; rm *.mod


