################################################################################
####################      Makefile for exaMD library     #######################
################################################################################

include make.inc

#####################
### Input sources ###
#####################

SOURCES:= \
	user_func.f90 exaMD.f90 

OBJECTS:=$(patsubst %.f90,%.o,$(SOURCES))

#########################
### Library locations ###
#########################

INCLUDES:=-ISRC/ -I$(METIS_DIR)/include
LIBS:=-LLIB/ -lexaMD -L$(METIS_DIR)/lib -lmetis 

###########################
### Compilation targets ###
###########################

default: exaMD

exaMD: lib $(OBJECTS)
	$(F90) $(FLAGS) -o $@ $(OBJECTS) $(LIBS)

%.o: %.f90 Makefile 
	$(F90) $(FLAGS) -c -o $@ $(INCLUDES) $<

lib: 
	mkdir -p LIB ; cd SRC ; $(MAKE) 

clean:
	rm -f LIB/* SRC/*.o SRC/*.mod *.o *.mod exaMD
