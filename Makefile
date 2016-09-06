# -*- Makefile -*- 
#
# Makefile for AnalyzeTool-x.x
# 
# by ACV
#
#


################ edit this section if you need ###########

BINDIR	:= $(HOME)/bin
EXEC	:= Analyze

# --- compilers and flags:
FC 	:= /scratch/cantal/opt/hdf5-1.8.17/hdf5/bin/h5fc
CC      := /scratch/cantal/opt/hdf5-1.8.17/hdf5/bin/h5cc
FC77    := gfortran
FFLAGS  := -lStARTLib -lcfitsio -ffree-line-length-none -static-libgfortran
LFLAGS  := $(FFLAGS) 
CFLAGS	:= -I $(HOME)/codes/lib -L $(HOME)/codes/lib
CCFLAGS := 
LINKER  := /scratch/cantal/opt/hdf5-1.8.17/hdf5/bin/h5fc
RM 	:= rm -f -v

##########################################################

# --- sources are here:
#include ./Makefile.sources

default: 
	$(LINKER) $(CFLAGS) $(LFLAGS) -c AnalysisTool.f90
	$(LINKER) $(CFLAGS) $(LFLAGS) AnalysisTool.o -o $(BINDIR)/$(EXEC)

clean:
	$(RM) *.o *.mod

tidy:	clean
	$(RM) $(BINDIR)/$(EXEC)

new:	clean default
 

 



