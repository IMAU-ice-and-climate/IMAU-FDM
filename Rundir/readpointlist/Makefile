FC = mpifort #gfortran
CC = cc

.SUFFIXES: .f90 .o

all: readpointlist.x distribute_points.x

readpointlist.x: readpointlist.o 
	$(FC) $(OFLAGS) -o $@ readpointlist.o 
	rm readpointlist.o 
	

distribute_points.x: sleep.o distribute_points.o 
	$(FC) $(OFLAGS) -o $@ sleep.o distribute_points.o 
	rm sleep.o distribute_points.o 
	
.f90.o:
	$(FC) $(FFLAGS) -c -o $*.o $<
.c.o:
	$(CC) -s -O2 -c -o $*.o $<
