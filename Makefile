PROCESSOR := $(shell uname -m)

F90=gfortran
FFLAGS=-g -C -O3 -ffree-form -fcheck=all -I/opt/local/include -fbounds-check
FFLAGS2=$(FFLAGS)
LDFLAGS=-L/opt/local/lib -lnetcdff -framework vecLib

.PHONY: clean 1d_test test

SUBDIR = hrefine/n5/

SOURCES= nDGsweep.f90 \
				 positivityLimit.f90

MODULES = testParameters.f90 \
					nDGmod.f90

SOURCES2= nDGmod.f90

OBJECTS=$(SOURCES:.f90=.o)
MODOBJ=$(MODULES:.f90=.o)

OBJECT2=$(SOURCES2:.f90=.o)

all: $(SOURCES) nodal_test.out

1d_test: nodal_test.out
	./nodal_test.out 2>&1 | tee screen.out
	cp screen.out _ndgunlim/$(SUBDIR)
	cp screen.out _ndgzhshu/$(SUBDIR)
	cp screen.out _matrunc/$(SUBDIR)

output: nodal_test.out
	./nodal_test.out

nodal_test.out: $(MODOBJ) $(OBJECTS) nodal_execute.f90
	$(F90) $(FFLAGS) $(MODOBJ) $(OBJECTS) nodal_execute.f90 -o $@ $(LDFLAGS)

quadNodes: quadNodes.out
	./quadNodes.out

quadNodes.out: $(OBJECT2) writeQuadData.f90
	$(F90) $(FFLAGS) $(OBJECT2) writeQuadData.f90 -o $@ $(LDFLAGS)

clean:
	rm -f *.o *.out *.mod *.nc

%.o : %.f90
	$(F90) -c $(FFLAGS) $<
