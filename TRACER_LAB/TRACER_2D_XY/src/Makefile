###################   Makefile   ####################
#                                                   #
# auteur : F.Voitus                                 #
# sujet  : Crée l'executable '../exe' avec './*'    #
#                                                   #
#####################################################

#-----------  Definition des variables  ------------#
COMPIL = mpif90
#COMPIL = gfortran
FLAGS  = -ffree-line-length-none \
	 -fno-protect-parens     \
	 -fcheck=bounds          \
	 -funroll-loops          \
	 -march=native           \
	 -fbackslash             \
	 -ffast-math             \
	 -fopenmp                \
	 -g3                     \
	 -O2                     \
         -mcmodel=large

OBJS   = modparam.o     \
         modfft.o                  \
	 moddiag.o               \
	 modsetup.o             \
	 modmpdta.o           \
	 modslag.o               \
         modview.o              \
	 modstep.o
EXE    = ../tracer_xy
#---------------------------------------------------#

#-----------  Creation de l'executable  ------------#
all: $(OBJS) main.o
	@$(COMPIL) $(FLAGS) -o $(EXE) $(OBJS) main.o -llapack -lblas -L. -lgft
	@echo -n "\n Succès de la compilation à " 
	@date +%R
	@echo " "
#---------------------------------------------------#

#------------  Compilation des modules  ------------#
main.o: main.F90 $(OBJS)
	@echo -n "main.o      : "
	@$(COMPIL) $(FLAGS) -c main.F90
	@echo "    ok"

modparam.o : modparam.f90
	@echo  "\n Édition de liens \n"
	@echo -n "modparam.o  : "
	@$(COMPIL) $(FLAGS) -c modparam.f90
	@echo "    ok"

modfft.o : modfft.f90 modparam.o
	@echo -n "modfft.o    : "
	@$(COMPIL) $(FLAGS) -c modfft.f90 -L. -lgft
	@echo "    ok"

moddiag.o : moddiag.f90 modparam.o
	@echo -n "moddiag.o   : "
	@$(COMPIL) $(FLAGS) -c moddiag.f90
	@echo "    ok"

modsetup.o : modsetup.f90 modfft.o modparam.o
	@echo -n "modsetup.o  : "
	@$(COMPIL) $(FLAGS) -c modsetup.f90
	@echo "    ok"

modslag.o : modslag.f90 modsetup.o modparam.o
	@echo -n "modslag.o   :"
	@$(COMPIL) $(FLAGS) -c modslag.f90
	@echo "     ok"

modmpdta.o : modmpdta.f90 modsetup.o modparam.o
	@echo -n "modmpdta.o  : "
	@$(COMPIL) $(FLAGS) -c modmpdta.f90
	@echo "    ok"

modview.o : modview.f90 moddiag.o modsetup.o modparam.o
	@echo -n "modview.o   : "
	@$(COMPIL) $(FLAGS) -c modview.f90
	@echo "    ok"

modstep.o : modstep.f90 modview.o moddiag.o modsetup.o modfft.o modparam.o
	@echo -n "modstep.o   : "
	@$(COMPIL) $(FLAGS) -c modstep.f90
	@echo "    ok"

#---------------------------------------------------#

#-------------------  Nettoyage  -------------------#
clean:
	@clear
	@clear
	@rm -f *~ *.o mod_*.mod 
	@rm -f $(EXE) ../*~
	@echo "Repertoire propre"
	@echo " "

proper:
	@clear
	@clear
	@rm -f *~ *.o *.mod *.a
	@rm -f ../$(EXE) ../*~ ../res/*
	@echo "Repertoire propre"
	@echo " "
#---------------------------------------------------#
