# Compiler 
CC = mpicxx
# Include
IFLAGS = /opt/intel/impi/include64/
#libraries
#LFLAGS = -shared-intel
LFLAGS = -floop-parallelize-all -ftree-parallelize-loops=2
all: proj
proj: proj.o 
	$(CC) $(LFLAGS) proj.o  -o proj -lm 
proj.o: proj.cpp
	$(CC) $(LFLAGS) -c  proj.cpp 
clean:
	find . -exec touch {} \;
	rm -rf *o proj
