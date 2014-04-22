PROCESSOR := $(shell uname -m)

F90=gfortran
FFLAGS=-g -C -O3 -ffree-form -fcheck=all -I/opt/local/include -fbounds-check
FFLAGS2=$(FFLAGS)
LDFLAGS=-L/opt/local/lib -lnetcdff -framework vecLib

.PHONY: clean nodal test

SOURCES= nDGmod.f90 nDGsweep.f90
SOURCES2= nDGmod.f90
OBJECTS=$(SOURCES:.f90=.o)
OBJECT2=$(SOURCES2:.f90=.o)

all: $(SOURCES) nodal_test.out

nodal: nodal_test.out
	./nodal_test.out

nodal_test.out: $(OBJECTS) nodal_execute.f90
	$(F90) $(FFLAGS) $(OBJECTS) nodal_execute.f90 -o $@ $(LDFLAGS) 

quadNodes: quadNodes.out
	./quadNodes.out

quadNodes.out: $(OBJECT2) writeQuadData.f90
	$(F90) $(FFLAGS) $(OBJECT2) writeQuadData.f90 -o $@ $(LDFLAGS)

clean:
	rm -f *.o *.out *.mod *.nc

%.o : %.f90
	$(F90) -c $(FFLAGS) $<

