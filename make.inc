################################################################################
####################   make.inc file for exaMD library   #######################
################################################################################

#######################
### Compiler set-up ###
#######################

# Compiler command
F90=mpif90

# Compiler flags
# The PRC macro sets the precision of the floating-point numbers and 
# the IPRC macro sets the precision of certain integers within the code. 
FLAGS= -O3 -cpp -DPRC=8 -DIPRC=4

# === Switches to change code ===

# Uncomment next line to print the calculated forces^2 for each block
#FLAGS+=-DPRINTALLFORCES

# Uncomment next line to use the diagonal connections 
# as well as the orthogonal ones in METIS partitioning
#FLAGS+=-DDODIAGS

#########################
### Library locations ###
#########################

# The location of the METIS library, used by the driver program.
#METIS_DIR:=$(HOME)/ParMETIS/Intel
METIS_DIR:=$(HOME)/ParMETIS/GNU
