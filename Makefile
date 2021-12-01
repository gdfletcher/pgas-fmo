#  PGAS-FMO proxy app Makefile 

 FC=mpifort
#FC=mpif90
 FLAGS=-O3 -w

 SRC=pgas-fmo.f90
#SRC=pgas-fmo.0.0.f90

all: pgas-fmo.x  

pgas-fmo.x: src/$(SRC)  
	$(FC) $(FLAGS) -o bin/pgas-fmo.x src/$(SRC) 

clean:
	rm *.mod


