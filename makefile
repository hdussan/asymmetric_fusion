#
#
#   makefile for asymmetric fusion project
#
#
OBJ = main.o potentials.o folded_potential.o semi_classical.o

CFLAGS= -c -O3 
COMPILE=g++
GSLFLAGS=-lgsl -lgslcblas
#
BoostPath=/usr/local/include/boost
InclBoost=-I$(BoostPath)
#
EIGENPATH=/usr/local/include/Eigen/
InclEigen=-I$(EIGENPATH)

#  List of header files for dependencies
MainHeader = a_fusion.h
PotHeader = potentials.h 
FoldedHeader = folded_potential.h
SemiClassicHeader = semi_classical.h
MinimalHeader = minimal.h
# Prepare for Dependencies
MainFile     = main.o
PotFiles     = potentials.o
FoldedFile   = folded_potential.o
SemiClassicFile = semi_classical.o

# Compilation
a_fusion:$(OBJ)
	$(COMPILE) $(GSLFLAGS) -o a_fusion $(OBJ)
$(OBJ):%.o: %.cpp
	$(COMPILE) $(CFLAGS) $(InclEigen) $(InclBoost) $< -o $@

# Dependencies
$(MainFile):  $(MainHeader)
$(PotFiles):  $(PotHeader)
$(FoldedFile):  $(FoldedHeader)
$(SemiClassicFile): $(SemiClassicHeader) 

#Standard cleaning
clean:
	rm *.o
cleanAll:
	rm a_fusion  *.o