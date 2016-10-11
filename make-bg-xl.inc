################################################################################
######   make.inc file for exaMD library using the Blue Gene XL compiler #######
################################################################################

#######################
### Compiler set-up ###
#######################

# Compiler command
F90=mpixlf90

# Compiler flags
# The PRC macro sets the precision of the floating-point numbers and 
# the IPRC macro sets the precision of certain integers within the code. 
FLAGS= -O3 -qsuffix=f=f90:cpp=f90 -qfree=f90 -WF,-DPRC=8 -WF,-DIPRC=4 -qstrict

# === Switches to change code ===

# Uncomment next line to print the calculated forces^2 for each block
#FLAGS+=-WF,-DPRINTALLFORCES

# Uncomment next line to use the diagonal connections 
# as well as the orthogonal ones in METIS partitioning
#FLAGS+=-WF,-DDODIAGS

#########################
### Library locations ###
#########################

# The location of the METIS library, used by the driver program.
METIS_DIR:=$(HOME)/ParMETIS/XL
