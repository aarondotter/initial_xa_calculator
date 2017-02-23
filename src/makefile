# This is the makefile for intial_xa_calculator

#MESA_DIR = 
#if you don't set this as an environment variable, set it here

#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

#LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lnet
LOAD_OTHER = -L$(MESA_LIB_DIR) -lnet  -leos -lrates -lchem $(LOAD_MESA_NUMERICS)

#################################################################

# STEP 3: build

CALC = initial_xa_calculator
CALC_OBJS = initial_xa_calculator.o

all : $(CALC)

$(CALC) : $(CALC_OBJS)
	$(LOADER) $(FCopenmp) -o ../$(CALC) $(CALC_OBJS) $(LOAD_OTHER)
	
#################################################################
   
initial_xa_calculator.o: initial_xa_calculator.f90
	$(WORK_COMPILE) $<

clean:
	rm -rf initial_xa_calculator initial_xa_calculator.o 
