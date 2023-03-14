
FC90 = gfortran
CC = gcc

OPTS90 = -O3

CCOPTS = -O3 -Wall -std=c99

LEGACY =


# Linkage rules

all : fdpic1

fdpic1 : fdpic1.o fdpush1.o dtimer.o dfield1.o edist.o
	$(FC90) $(OPTS90) -o fdpic1 fdpic1.o fdpush1.o dtimer.o dfield1.o edist.o

# Compilation rules

dfield1.o : dfield1.f
	$(FC90) $(OPTS90) -o dfield1.o -c dfield1.f

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fdpush1.o : dpush1.f
	$(FC90) $(OPTS90) -o fdpush1.o -c dpush1.f


edist.o: edist.f90
	$(FC90) $(OPTS90) -o edist.o -c edist.f90

fdpic1.o : dpic1.f90 
	$(FC90) $(OPTS90) -o fdpic1.o -c dpic1.f90


clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fdpic1
