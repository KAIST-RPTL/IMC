#SHELL = /bin/bash

COMPILERC   = icc
COMPILERF   = mpiifort
FLAG  		= -fopenmp -mkl -Mmodules ./obj/ -O2
FLAGC       = -fopenmp -mkl -O2
MPIFLAG	  	= -mcmodel=medium -heap-arrays -mkl

CODE 		= IMC.out 

SRC_PATH := ./
OBJ_PATH := ./obj/

include $(SRC_PATH)/FILE.list

FOBJS :=$(patsubst %,$(OBJ_PATH)%,$(FOBJS))
COBJS :=$(patsubst %,$(OBJ_PATH)%,$(COBJS))

.SUFFIXES: .f90 .o .c
$(OBJ_PATH)%.o : $(SRC_PATH)%.c
				$(COMPILERC) $(FLAGC) -o $@ -c $<
$(OBJ_PATH)%.o : $(SRC_PATH)%.f90
				$(COMPILERF) $(FLAG) $(BUG) $(PROFILE) -o $@ -c $<
$(CODE) : $(FOBJS) $(COBJS)
	$(COMPILERF) $(FLAG) $(PROFILE) -o  $(CODE) $(FOBJS) $(COBJS)

# Make clean
clean :
	@echo "Cleaning files"
	@rm  $(OBJ_PATH)* -f
	@rm  *.mod -f
	@rm  $(CODE) -f
