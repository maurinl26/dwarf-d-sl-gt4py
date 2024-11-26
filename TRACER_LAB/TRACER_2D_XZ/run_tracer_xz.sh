
#!/bin/bash
#####################  SHELL LANÇANT ET CLASSANT LES EXPÉRIENCES  #####################
#                                                                                     #
# auteurs : L.Auger, F.Voitus, C.Colavolpe                                            #
# sujet   : Lance le modèle après avoir compilé et classe les resultats               #
#                                                                                     #
#######################################################################################

#---------------------  Definition des variables  ---------------------#
EXP="TEST"

MACHINE="lxgmap35:"
CHEMIN=$PWD
#CHEMIN=${MACHINE}/home/voitus/Desktop/NEW_DYN/TIME/SI/MODEL/2D_SLICI_BIGWADV

SCH=$(grep "SYSTEM=" init)
SCH=${SCH:12:2}
#EXP=${EXP}"_"${SCH}
#---------------------------------------------------------------------#

#----------------  Compilation et nettoyage en local  ----------------#
cd src/
touch main.F90
make
cd ..

if [ -d exp/${EXP} ];then
    rm -rf  exp/${EXP}
fi
mkdir -p exp/${EXP}/res
cp init ./exp/${EXP}
#---------------------------------------------------------------------#

#----------------------------  Execution  ----------------------------#
if [ -d res ];then
    rm -rf  res
fi
mkdir res
./tracer_xz
mv ./res/* exp/${EXP}/res
rm -rf res
#---------------------------------------------------------------------#


