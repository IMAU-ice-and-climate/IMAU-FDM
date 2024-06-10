FC = ftn
CC = cc

.SUFFIXES: .f90 .o

# intel
#OFLAGS  = -g
#FFLAGS  = -g -c -ffixed-line-length-132 -w -s -O -fbounds-check
# cray
OFLAGS = 
FFLAGS = -c -N 132 -O2 -h bounds

#OBJECTS = \

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
