##################################################################
#                        Makefile for f90                         #
###################################################################

########################     compiler     #########################

F90 = mpiifx
F77 = mpiifx
#F90 = gfortran
#F77 = gfortran

########################  compiler flags  #########################

F90FLAGS = -c -O0 -g -check bounds -traceback
#F90FLAGS = -c -warn -CB -debug extended
F77FLAGS = -c
LFLAGS   =
#PREP = scalasca -inst
#PREP = skin

########################  objects alpha   #########################

INIT    = .
SRCDIR  = $(INIT)
OBJ     = $(INIT)
OBJDIR  = $(OBJ)
CALCDIR = $(INIT)

OBJECTS = $(OBJ)/declaration.o\
          $(OBJ)/unify_map.o

########################      build       #########################

unify_map : printmsgA $(OBJECTS)
	@echo
	@echo Linking...
	$(PREP) $(F90) -o $@ $(OBJECTS) $(LFLAGS)
	@echo
	@echo unify_map built, congratulations.
	@echo

########################     compile      #########################

$(OBJDIR)/declaration.o : $(SRCDIR)/declaration.f90 $(SRCDIR)/Makefile
	@echo compiling declaration.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/declaration.f90

$(OBJDIR)/unify_map.o : $(SRCDIR)/unify_map.f90 $(SRCDIR)/Makefile
	@echo compiling unify_map.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/unify_map.f90

########################      message     #########################

printmsgA :
	@echo
	@echo Building ...
	@echo F77 Compiler flags : $(F77FLAGS)
	@echo F90 Compiler flags : $(F90FLAGS)
	@echo Linker   flags : $(LFLAGS)
	@echo Prepend  flags : $(PREP)

########################      clean       #########################

clean:
	@find . \( -name '*.o' -o -name '*.mod' -o -name 'unify_map' \) -exec rm -f {} \;

########################   end of file    #########################