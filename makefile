#
#
#   makefile for asymmetric fusion project
#
#
OBJ = main.o potentials.o folded_potential.o semi_classical.o

CFLAGS= -c -O3 
COMPILE=g++
THREADFLAG=-lpthread
#
BoostPath=/usr/local/include/boost
InclBoost=-I$(BoostPath)
# 
#   boost threads
#
IncludePath=/usr/local/include
InclInclude=-I$(IncludePath)
LibPath=/usr/local/lib
LibLink=-L$(LibPath)
BoostThreadFLAG=-lboost_system -lboost_thread-mt

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
	$(COMPILE) $(InclInclude) $(LibLink) $(BoostThreadFLAG) -o a_fusion $(OBJ)
$(OBJ):%.o: %.cpp
	$(COMPILE) $(CFLAGS) $(InclBoost) $< -o $@

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