#!/bin/bash
###################   ee_3D_eta   ###################
#                                                   #
# auteur : C.Colavolpe                              #
# sujet  : - Compile via './src/*'                  #
#          - Execute './exe'                        #
#          - Classe des sorties dans './xp/*'       #
#                                                   #
#####################################################

#-----------  Definition des variables  ------------#
XP=TEST

NP=4
NC=1
#---------------------------------------------------#

#------------------  Compilation  ------------------#
rm -rf "a.out"

cd src/
echo -e "\n Mise à jour du programme \n"
touch modparam.f90
make
cd ..

if [ ! -f "slici_3d" ];then
    echo -e "\n Erreur de compilation \n "
    exit
fi
#---------------------------------------------------#

#-------  Nettoyage précédentes expériences  -------#
if [ -d exp/${XP} ];then
#    echo -e "supprimer le dossier '"exp/${XP}"' ? (o/n) "
#    read ON1
#    if [ ${ON1} == 'n' ];then
#	exit 
#    fi
    rm -rf  exp/${XP}
fi

mkdir -p exp/${XP}/res/XY
mkdir -p exp/${XP}/res/XZ
cp init ./exp/${XP}
cp ./src/modparam.f90 ./exp/${XP}

if [ -d res ];then
     echo -e "supprimer le dossier res ? (o/n) "
     read ON2
     if [ ${ON2} == 'n' ];then
	 exit 
     fi
     rm -rf res${NB}
fi
#---------------------------------------------------#

#-------------------  Execution  -------------------#
echo -e " \n Expérience      : "${XP}
echo -e " dans le dossier : "res${NB}
mkdir -p res/XY
mkdir -p res/XZ
export OMP_NUM_THREADS=${NC}
ulimit -s unlimited
mpirun -q -np ${NP} tracer_xy
mv ./res/* exp/${XP}/res/
rm -rf res
echo -e " Fin de l'expérience : "${XP}"\n"
#---------------------------------------------------#


