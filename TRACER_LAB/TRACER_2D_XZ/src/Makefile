#######################  MAKEFILE ########################
# auteurs : F.Voitus                                     #
##########################################################

#---------------------  Definition des varialbes  ---------------------#
COMPIL = gfortran
FLAGS  = -ffree-line-length-none \
	-fno-protect-parens      \
	-g3 -fcheck=bounds       \
	-funroll-loops           \
	-march=native            \
	-fbackslash              \
	-ffast-math -O3          \
	-flto                    
#        -fno-strict-aliasing
EXEDIR = ..
OBJS   = modshare.o       	 \
	 modfft.o         	 \
	 modgrad.o               \
	 modsetup.o              \
	 moddiag.o               \
	 modslag.o               \
	 modmpdta.o              \
         modview.o               \
	 modstate.o
EXE    = tracer_xz
#---------------------------------------------------------------------#

#-------------  Compilation et creation de l'executable  -------------#

all: $(OBJS) main.o
	@$(COMPIL) $(FLAGS) -o $(EXEDIR)/$(EXE) $(OBJS) main.o -L. -lgft -llapack
	@echo -n "\n Succès de la compilation à "
	@date +%R
	@echo " "

#---------------------------------------------------#
#------------  Compilation des modules  ------------#

main.o: main.F90 $(OBJS) 
	@echo -n "main.o : "
	@$(COMPIL) $(FLAGS) -c main.F90
	@echo "         ok"

modshare.o : modshare.f90
	@echo  "\n Début de la compilation \n"
	@echo -n "modshare.o : "
	@$(COMPIL) $(FLAGS) -c modshare.f90
	@echo "     ok"
modfft.o : modfft.f90 modshare.o
	@echo -n "modfft.o : "
	@$(COMPIL) $(FLAGS) -c modfft.f90 -L. -lgft
	@echo "       ok"
modgrad.o : modgrad.f90 modfft.o modshare.o
	@echo -n "modgrad.o : "
	@$(COMPIL) $(FLAGS) -c modgrad.f90
	@echo "      ok"
moddiag.o : moddiag.f90 modshare.o
	@echo -n "moddiag.o : "
	@$(COMPIL) $(FLAGS) -c moddiag.f90
	@echo "      ok"
modsetup.o : modsetup.f90 moddiag.o modgrad.o modshare.o
	@echo -n "modsetup.o : "
	@$(COMPIL) $(FLAGS) -c modsetup.f90
	@echo "     ok"
modslag.o : modslag.f90 modsetup.o modshare.o
	@echo -n "modslag.o : "
	@$(COMPIL) $(FLAGS) -c modslag.f90
	@echo "      ok"
modmpdta.o : modmpdta.f90 modsetup.o modshare.o
	@echo -n "modmpdta.o : "
	@$(COMPIL) $(FLAGS) -c modmpdta.f90
	@echo "     ok"	
modview.o : modview.f90 moddiag.o modsetup.o modshare.o
	@echo -n "modview.o : "
	@$(COMPIL) $(FLAGS) -c modview.f90
	@echo "      ok"
modstate.o : modstate.f90 moddiag.o modsetup.o modgrad.o modshare.o
	@echo -n "modstate.o : "
	@$(COMPIL) $(FLAGS) -c modstate.f90
	@echo "     ok"

#---------------------------------------------------------------------#
#--------------------  Nettoyage des repertoires  --------------------#
clean:
	@clear
	@clear
	@rm -f *~ *.o mod_*.mod
	@rm -f ../$(EXE) ../*~
	@echo "Repertoire propre"
	@echo " "
#---------------------------------------------------------------------#
#--------------------  Nettoyage des repertoires  --------------------#
proper:
	@clear
	@clear
	@rm -f *~ *.o *.mod *.a
	@rm -f ../$(EXE) ../*~ ../res/*
	@echo "Repertoire propre"
	@echo " "
#--------------------------------------------------------------------#
