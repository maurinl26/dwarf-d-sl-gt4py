###################   Makefile   ####################
#                                                   #
# auteur : F.Voitus                                 #
# sujet  : Crée l'executable '../exe' avec './*'    #
#                                                   #
#####################################################

#-----------  Definition des variables  ------------#
COMPIL = mpif90
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

OBJS   = modparam.o              \
         modfft.o                \
	 moddiag.o               \
	 modstruc.o              \
	 modcomm.o               \
	 modmpdta.o
	 modslag.o               \
         modview.o               \
	 moddt.o
EXE    = ../tracer_3d
#---------------------------------------------------#

#-----------  Creation de l'executable  ------------#
all: $(OBJS) main.o
	@$(COMPIL) $(FLAGS) -o $(EXE) $(OBJS) main.o -llapack -lblas -L. -lgft
	@echo -en "\n Succès de la compilation à " 
	@date +%R
	@echo " "
#---------------------------------------------------#

#------------  Compilation des modules  ------------#
main.o: main.F90 $(OBJS)
	@echo -n "main.o      : "
	@$(COMPIL) $(FLAGS) -c main.F90
	@echo "    ok"
	@echo -e "\n Début de la compilation"

modparam.o : modparam.f90
	@echo -e "\n Édition de liens \n"
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

modstruc.o : modstruc.f90 modlinalg.o modfft.o modparam.o
	@echo -n "modstruc.o  : "
	@$(COMPIL) $(FLAGS) -c modstruc.f90
	@echo "    ok"

modcomm.o : modcomm.f90 modstruc.o modparam.o
	@echo -n "modcomm.o   : "
	@$(COMPIL) $(FLAGS) -c modcomm.f90
	@echo "    ok"

modslag.o : modslag.f90 modcomm.o modstruc.o modparam.o
	@echo -n "modslag.o   :"
	@$(COMPIL) $(FLAGS) -c modslag.f90
	@echo "     ok"

modmpdta.o : modmpdta.f90 modcomm.o modstruc.o modparam.o
	@echo -n "modmpdta.o : "
	@$(COMPIL) $(FLAGS) -c modmpdta.f90
	@echo "    ok"

modview.o : modview.f90 moddiag.o modstruc.o modparam.o
	@echo -n "modview.o   : "
	@$(COMPIL) $(FLAGS) -c modview.f90
	@echo "    ok"

moddt.o : moddt.f90 modview.o moddiag.o modcomm.o modstruc.o modfft.o modparam.o
	@echo -n "moddt.o     : "
	@$(COMPIL) $(FLAGS) -c moddt.f90
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
