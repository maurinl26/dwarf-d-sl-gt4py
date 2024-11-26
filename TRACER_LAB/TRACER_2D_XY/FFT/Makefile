####################   MAKEFILE COMPILANT LES FFT DU MODELE JOUET  #################### 
#                                                                                     #
# Auteurs : J.Chergui, L.Auger, F.Voitus, C.Colavolpe                                 #
# sujet   : Compile l'ensemble des modules FFT pour construire un modele jouet, et    #
#           mets à jour dans le dossier .src/                                         #
#######################################################################################

#---------------------  Definition des varialbes  ---------------------#
COMPIL  = gfortran
FLAGS  = -ffree-line-length-none \
	 -fno-protect-parens     \
	 -fcheck=bounds          \
	 -funroll-loops          \
	 -march=native           \
	 -fbackslash             \
	 -ffast-math             \
	 -g3                     \
	 -O2                     \
         -mcmodel=large
AR      = ar
ARFLAGS = rv
GFTLIB  = -L$(GFTLIBPATH) -lgft
LIBS    = $(GFTLIB)
LIBNAME = gft
LIB	= lib$(LIBNAME).a
OBJS	= $(LIB)(GFT_common.o)  \
          $(LIB)(GFT_jmfft.o)   \
          $(LIB)(GFT_set.o)     \
          $(LIB)(GFT_do.o)      \
          $(LIB)(GFT_end.o)     \
          $(LIB)(GFT.o)
#---------------------------------------------------------------------#

.SUFFIXES : .f .f90 .o .a

#-------------------  Compilation des modules FFT  -------------------#
all: $(OBJS)
	@rm -f ../src/GFT*
	@mv libgft.a ../src/
	@mv *.mod ../src/
#---------------------------------------------------------------------#

#-------------------  Compilation des modules FFT  -------------------#
clean:
	 rm -f *.o *gft.a *.mod *~
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#

.f.a:
	@$(COMPIL) -c $(FLAGS1) $(FLAGS2) $(*).f
	@$(AR) $(ARFLAGS) $(LIB) $(*).o
	@rm -f $(*).o
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
.f90.a:
	@$(COMPIL) -c $(FLAGS) $(*).f90
	@$(AR) $(ARFLAGS) $(LIB) $(*).o 
	@rm -f $(*).o
#---------------------------------------------------------------------#
