# Makefile for PROTEUS

#--- Project -------------------------------------------

PROJECT = Jupiter/Nakamura_22_Case1

#--- Compiler ------------------------------------------

FC      = /usr/local/bin/gfortran
FFLAGS  = -O3

#--- Source --------------------------------------------

SRCS   = src/e__main.f90 \
         src/v__tdec.f90 \
         src/c__prm.f90 \
         src/p__io.f90 \
         src/p__search.f90 \
         src/v__Venus.f90 \
         src/v__Earth.f90 \
         src/v__Mars.f90 \
         src/v__Jupiter.f90 \
         src/p__EUVAC.f90 \
         src/p__EUV.f90 \
         src/p__UV.f90 \
         src/p__eddy_diffusion.f90 \
         src/p__molecular_diffusion.f90 \
         src/p__photolysis_rate.f90 \
         src/p__photoionization_rate.f90 \
         src/p__photochem_rate.f90 \
         src/p__photochem_transport.f90 \
         src/p__vertical_diffusion.f90 \
         src/p__chemical_scheme.f90 \
         src/p__photochem_scheme.f90 \
         $(PROJECT)/p__PROTEUS.f90 \
         $(PROJECT)/v__in.f90 

SRCSF  = src/HEUVAC.F

#--- Target --------------------------------------------

TARGET  = ../PROTEUS_exe
OBJS    = $(SRCS:.f90=.o)
OBJSF   = $(SRCSF:.F=.o)
INCDIR  = src
INCDIR2 = $(PROJECT)

#--- Rules ---------------------------------------------

.PHONY: all clean
.SUFFIXES : .o .f90 .F

#--- Compile -------------------------------------------

all: MSG1 $(INCDIR) $(INCDIR2) $(TARGET) MSG2
MSG1:
	@echo Compiling PROTEUS...
MSG2:
	@echo Finished compiling PROTEUS.
	@echo Executable file "$(TARGET)" was successfully generated.
	@echo To run PROTEUS, go to \"..\" directory and type \"./PROTEUS_exe\".
$(TARGET) : $(OBJS) $(OBJSF)
	$(FC) $(FFLAGS) -o $@ $^
$(INCDIR):
	@mkdir -p $(INCDIR)
$(INCDIR2):
	@mkdir -p $(INCDIR2)
$(OBJS): %.o: ../%.f90
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJSF): %.o: ../%.F
	$(FC) $(FFLAGS) -c $< -o $@

#--- Dependencies --------------------------------------

src/e__main.o: src/v__tdec.o \
               src/c__prm.o \
               src/p__io.o \
               src/p__search.o \
               src/v__Venus.o \
               src/v__Earth.o \
               src/v__Mars.o \
               src/v__Jupiter.o \
               src/p__EUVAC.o \
               src/p__EUV.o \
               src/HEUVAC.o \
               src/p__UV.o \
               src/p__eddy_diffusion.o \
               src/p__molecular_diffusion.o \
               src/p__photolysis_rate.o \
               src/p__photoionization_rate.o \
               src/p__photochem_rate.o \
               src/p__photochem_transport.o \
               src/p__vertical_diffusion.o \
               src/p__chemical_scheme.o \
               src/p__photochem_scheme.o \
               $(PROJECT)/p__PROTEUS.o \
               $(PROJECT)/v__in.o
src/v__tdec.o:
$(PROJECT)/v__in.o: src/v__tdec.o src/p__search.o
src/v__Earth.o: src/v__tdec.o src/p__search.o
src/v__Venus.o: src/v__tdec.o src/p__search.o
src/v__Mars.o: src/v__tdec.o src/p__search.o
src/v__Jupiter.o: src/v__tdec.o src/p__search.o
src/c__prm.o: src/v__tdec.o 
src/p__search.o: src/v__tdec.o
$(PROJECT)/p__PROTEUS.o:
src/p__EUVAC.o: src/v__tdec.o src/p__search.o src/c__prm.o
src/p__EUV.o: src/v__tdec.o src/p__search.o src/c__prm.o src/HEUVAC.o src/p__photoionization_rate.o
src/p__UV.o: src/v__tdec.o src/p__search.o src/c__prm.o src/p__photolysis_rate.o $(PROJECT)/p__PROTEUS.o
src/p__eddy_diffusion.o: src/v__tdec.o src/p__search.o
src/p__molecular_diffusion.o: src/v__tdec.o src/p__search.o
src/p__vertical_diffusion.o:
src/p__photochem_rate.o: src/v__tdec.o $(PROJECT)/p__PROTEUS.o src/p__search.o src/p__photolysis_rate.o src/p__photoionization_rate.o
src/p__photochem_transport.o: src/v__tdec.o src/p__search.o src/p__eddy_diffusion.o src/p__molecular_diffusion.o src/p__vertical_diffusion.o
src/p__chemical_scheme.o:
src/p__photochem_scheme.o: src/v__tdec.o src/p__search.o src/p__io.o src/p__chemical_scheme.o
src/HEUVAC.o:
src/p__io.o: src/v__tdec.o src/p__search.o

#--- Clean ---------------------------------------------

clean:
	rm -f $(OBJS) $(OBJSF) *.mod $(TARGET) ../src/*.mod

