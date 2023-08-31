#SHELL = /bin/bash

COMPILERC   = icc
COMPILERF   = mpiifort
FLAG  		= -march=core-avx2 -qopenmp -mcmodel=medium -heap-arrays -mkl -module ./obj/ -O2
FLAGC       = -march=core-avx2 -qopenmp -mcmodel=medium -heap-arrays -mkl -O2
#BUG			= -check all -traceback -fp-stack-check
MPIFLAG	  	= -mcmodel=medium -heap-arrays -mkl

### TAU UTIL ### 
#  #export PATH=$PATH:/home/guest/HyeonTae/Traiing/tau-2.28.1/bin
#  export TAU_PROFILE=/home/guest/HyeonTae/Training/tau-2.28.1/x86_64/lib/Makefile.tau-icpc-mpi-pdt-openmp
#  #TAU_PROFILE ?=/home/guest/HyeonTae/Training/tau-2.28.1/include/Makefile
#  include $(TAU_PROFILE) 
#  COMPILERF = TAU_PROFILE=$(TAU_PROFILE) /home/guest/HyeonTae/Training/tau-2.28.1/x86_64/bin/tau_f90.sh
#  COMPILERC = TAU_PROFILE=$(TAU_PROFILE) /home/guest/HyeonTae/Training/tau-2.28.1/x86_64/bin/tau_cc.sh
CODE 		= IMC.out 

#SRC_PATH := ../01_iMC_src/src_merge2/
#SRC_PATH := ../91_iMC_IFP/src_merge2/
SRC_PATH := ../source/DEPLETION/
#SRC_PATH := ../source/SOURCE_IFP/
OBJ_PATH := ./obj/

include $(SRC_PATH)/FILE.list
#include FILE.list

FOBJS :=$(patsubst %,$(OBJ_PATH)%,$(FOBJS))
COBJS :=$(patsubst %,$(OBJ_PATH)%,$(COBJS))

.SUFFIXES: .f90 .o .c
#.c.o:
#	$(COMPILERC) $(FLAG) $(BUG) $(PROFILE) -c  $*.c
#
#.f90.o:
#	$(COMPILERF) $(FLAG) $(BUG) $(PROFILE) -c  $*.f90
$(OBJ_PATH)%.o : $(SRC_PATH)%.c
				$(COMPILERC) $(FLAGC) -o $@ -c $<
$(OBJ_PATH)%.o : $(SRC_PATH)%.f90
				$(COMPILERF) $(FLAG) $(BUG) $(PROFILE) -o $@ -c $<
	
#%.o: %.mod

#a.f90.mod:
#	$(COMPILERF) $(FLAG) -c $<

# Make CODE
# all : $(CODE)
$(CODE) : $(FOBJS) $(COBJS)
	$(COMPILERF) $(FLAG) $(PROFILE) -o  $(CODE) $(FOBJS) $(COBJS)
	
# Make CODE - openmp
#omp : $(CODE)
#
#$(CODE) : $(FOBJS) $(COBJS)
#	$(COMPILERF) $(FLAG)  -o $(CODE) $(FOBJS) $(COBJS)
#
#noomp : $(CODE) 
#	$(COMPILERF) -o $(CODE) $(FOBJS) $(COBJS)

# Make clean
clean :
	@echo "Cleaning files"
	@rm  $(OBJ_PATH)* -f
	@rm  *.mod -f
	@rm  $(CODE) -f
