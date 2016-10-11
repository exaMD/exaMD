################################################################################
####################      Makefile for exaMD library     #######################
################################################################################

include ../make.inc

#####################
### Input sources ###
#####################

SOURCES:= \
	exaMD_constants_mod.f90 exaMD_module.f90 \
	exaMD_metis_interface.f90 exaMD_partitioning.f90 \
	exaMD_comms.f90 exaMD_comms_overlap.f90 exaMD_forces.f90 \
	exaMD_input_formats.f90 exaMD_utils.f90  

OBJECTS:=$(patsubst %.f90,%.o,$(SOURCES))

###########################
### Compilation targets ###
###########################

default: lib

lib: $(OBJECTS)
	ar -r ../LIB/libexaMD.a $(OBJECTS)

%.o: %.f90 Makefile 
	$(F90) $(FLAGS) -c -o $@ $<

clean:
	rm -f ../LIB/* *.o *.mod
