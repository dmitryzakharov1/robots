# Compiler 
CC = g++
# Include
#IFLAGS = /opt/intel/impi/include64/
#libraries
#LFLAGS = -shared-intel
#LFLAGS = -floop-parallelize-all -ftree-parallelize-loops=4 -fopenmp
LFLAGS =  -fopenmp
all: proj_omp
proj: proj_omp.o 
	$(CC) $(LFLAGS) proj_omp.o  -o proj_omp -lm 
proj_omp.o: proj_omp.cpp
	$(CC) $(LFLAGS) -c  proj_omp.cpp 
clean:
	find . -exec touch {} \;
	rm -rf *o proj_omp
