#!/bin/bash
###################   ee_3D_eta   ###################
#                                                   #
# auteur : C.Colavolpe                              #
# sujet  : Visalisation des résultats               #
#            ./resultat.sh FICHIER IT               #
#####################################################

#-----------  Definition des variables  ------------#
XP=HOT_EX
FQ=1
#---------------------------------------------------#

#-------------  Création des graphes  --------------#
IT=$(($2*${FQ}))
ncl -Q 'xp="'${XP}'"' 'plan="X"' 'step='${IT} $1
#ncl -Q 'xp="'${XP}'"' 'plan="Y"' 'step='${IT} $1
ncl -Q 'xp="'${XP}'"' 'plan="Z"' 'step='${IT} $1
#---------------------------------------------------#

#------------  Visionage et nettoyage  -------------#
if [ -f  T${IT}_Z.epsi ];then
  epstopdf T${IT}_X.epsi
#  epstopdf T${IT}_Y.epsi
  epstopdf T${IT}_Z.epsi
#  evince T${IT}_X.pdf T${IT}_Z.pdf T${IT}_Y.pdf  2> /dev/null
  evince T${IT}_X.pdf T${IT}_Z.pdf 2> /dev/null
fi

ITMAX=$(ls -A1 ${XP}/res/X/T* | wc -l)
echo -e "\n --> "$((${ITMAX}-1))" fichiers (hors état inital) dans le dossier "${XP}/res"\n"

rm -f *.epsi *.pdf
#---------------------------------------------------#
