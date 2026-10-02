#default fortran compiler
gnu_fc = gfortran -ffree-line-length-none
fc_legacy = gfortran -std=legacy -fallow-argument-mismatch -Ofast -ffree-line-length-none

#executable name
exec = test
volcano_driver = $(wildcard test_volcano.f90)

#flags
switchOPT = -O3 -flto -funroll-loops
switchDBG = -O0 -check all -warn all -fpe0 -u -traceback -warn nounused -g
switchPRO = $(switchOPT) -pg -traceback -g
switchOMP = $(switchOPT) -openmp

switch = $(switchOPT)

#objects that need fortran, mind the sequence
core_objs = opkda1.o opkda2.o opkdmain.o patmo_commons.o patmo_constants.o \
    	patmo_parameters.o patmo_utils.o patmo_rates.o patmo_reverseRates.o \
    	patmo_photo.o patmo_photoRates.o patmo_volc.o patmo_sparsity.o patmo_jacobian.o \
    	patmo_ode.o patmo.o
gnu_objs = $(core_objs) test.o
volcano_objs = $(core_objs) test_volcano.o

# default target
all: $(exec) $(if $(volcano_driver),test_volcano)

$(exec): $(gnu_objs)
	$(gnu_fc) $(switch) $(gnu_objs) -o $(exec)

test_volcano: $(volcano_objs)
	$(gnu_fc) $(switch) $(volcano_objs) -o test_volcano

.PHONY: check_volcano_transport
check_volcano_transport: $(core_objs) check_transport.f90
	$(gnu_fc) -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow $(filter-out patmo_volc.o,$(core_objs)) patmo_volc.f90 check_transport.f90 -o check_volcano_transport
	./check_volcano_transport

.PHONY: check_volcano_radiation
check_volcano_radiation: $(core_objs) check_radiation.f90
	$(gnu_fc) -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow $(core_objs) check_radiation.f90 -o check_volcano_radiation
	./check_volcano_radiation

# Rules for Fortran 77 files with legacy support
opkda1.o: opkda1.f
	$(fc_legacy) -c opkda1.f -o opkda1.o

opkda2.o: opkda2.f
	$(fc_legacy) -c opkda2.f -o opkda2.o

opkdmain.o: opkdmain.f
	$(fc_legacy) -c opkdmain.f -o opkdmain.o

debug: switch = $(switchDBG)
debug: nowarn = -nowarn
debug: all

profile: switch = $(switchPRO)
profile: nowarn = -nowarn
profile: all

omp: openmp
openmp: switch = $(switchOMP)
openmp: all

gfortran: fc = gfortran
gfortran: switch = -ffree-line-length-none
gfortran: all

gfortran_dbg: fc = gfortran
gfortran_dbg: switch = -fbacktrace -g
gfortran_dbg: switch += -ffpe-trap=zero,overflow,invalid
gfortran_dbg: switch += -fbounds-check -ffree-line-length-none
gfortran_dbg: all

coverage: fc = gfortran
coverage: switch = -ffree-line-length-none
coverage: switch += -fprofile-arcs -ftest-coverage
coverage: all

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec) test_volcano

test.o: test.f90 patmo.o patmo_photoRates.o

# Rebuild consumers when generated dimensions or physical constants change.
patmo_parameters.o: patmo_commons.o
$(filter-out opkda1.o opkda2.o opkdmain.o patmo_commons.o patmo_constants.o patmo_parameters.o,$(core_objs)): patmo_commons.o patmo_parameters.o patmo_constants.o
patmo_jacobian.o: patmo_utils.o

patmo_photoRates.o: patmo_parameters.o patmo_commons.o patmo_constants.o

patmo_volc.o: patmo_parameters.o patmo_utils.o patmo_rates.o patmo_reverseRates.o patmo_photoRates.o
patmo_ode.o: patmo_volc.o
patmo.o: $(filter-out patmo.o opkda1.o opkda2.o opkdmain.o,$(core_objs))

test_volcano.o: test_volcano.f90 patmo.o patmo_photo.o patmo_volc.o

#rule for f
%.o: %.f90
	$(gnu_fc) $(switch) -c $< -o $@

# Compile Fortran 77 files (legacy support)
%.o: %.f
	$(fc_legacy) -c $< -o $@
